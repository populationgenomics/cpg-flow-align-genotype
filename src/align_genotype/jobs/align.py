"""
FASTQ/BAM/CRAM -> CRAM: create Hail Batch jobs for (re-)alignment.
"""

import os.path
from textwrap import dedent
from typing import cast

from cpg_flow import filetypes, targets
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import BashJob

DRAGMAP_INDEX_FILES = ['hash_table.cfg.bin', 'hash_table.cmp', 'reference.bin']
PIPEFAIL = 'set -euo pipefail'


def align(
    sequencing_group: targets.SequencingGroup,
    job_attrs: dict,
    output_path: Path,
    markdup_metrics_path: str,
) -> list[BashJob]:
    """
    - if the input is 1 fastq pair, submits one alignment job.

    - if the input is a set of fastq pairs, submits multiple jobs per each pair,
    then submits a separate merge job.

    - if the input is a cram/bam

        - if number_of_shards_for_realignment is > 1, use bazam to shard inputs
        and align in parallel, then merge result together.

        - for dragmap, submit an extra job to extract a pair of fastqs from the cram/bam,
        because dragmap can't read streamed files from bazam.

    - deduplication (dupblaster) and CRAM conversion run in-stream within the
    alignment/merge job, so no separate markdup job is submitted.

    - nthreads can be set for smaller test runs on toy instance, so the job
    doesn't take entire 32-cpu/64-threaded instance.
    """

    batch_instance = hail_batch.get_batch()

    alignment_input = sequencing_group.alignment_input

    # consider setting the shard count to at least 2 when the input is large,
    # or else dragen-os can run out of memory and cause the VM to continually restart
    realignment_shards = config.config_retrieve(['workflow', f'{sequencing_group.sequencing_type}_realignment_shards'])

    sharded_fq = isinstance(alignment_input, filetypes.FastqPairs) and len(alignment_input) > 1
    sharded_bazam = (
        isinstance(alignment_input, filetypes.CramPath | filetypes.BamPath)
        and alignment_input.index_path
        and realignment_shards > 1
    )
    sharded = sharded_fq or sharded_bazam

    jobs = []
    sharded_align_jobs = []
    sorted_bams = []

    if not sharded:  # Just running one alignment job
        if isinstance(alignment_input, filetypes.FastqPairs):
            alignment_input = alignment_input[0]
        align_j, align_cmd, fifo_epilogue = _align_one(
            alignment_input=alignment_input,
            sequencing_group_name=sequencing_group.id,
            job_attrs=job_attrs,
            name_sort_for_merge=False,
        )
        jobs.append(align_j)
        merge_or_align_j = align_j

    # Aligning in parallel and merging afterwards
    else:
        fifo_epilogue = ''  # the merged command below never uses a fifo
        if sharded_fq:  # Aligning each lane separately, merging after
            # running alignment for each fastq pair in parallel
            fastq_pairs = cast('filetypes.FastqPairs', alignment_input)
            for pair in fastq_pairs:
                # dragmap command, name-sorting each lane (dedup happens after the merge):
                j, cmd, shard_fifo_epilogue = _align_one(
                    alignment_input=pair,
                    sequencing_group_name=sequencing_group.id,
                    job_attrs=job_attrs,
                    name_sort_for_merge=True,
                )
                if shard_fifo_epilogue:
                    cmd += f'\n{shard_fifo_epilogue}'
                j.command(cmd)
                sorted_bams.append(j.sorted_bam)
                sharded_align_jobs.append(j)

        elif sharded_bazam:  # Using BAZAM to shard CRAM
            for shard_number in range(realignment_shards):
                shard_string = f'{shard_number + 1}/{realignment_shards}'
                j, cmd, shard_fifo_epilogue = _align_one(
                    alignment_input=alignment_input,
                    sequencing_group_name=sequencing_group.id,
                    job_attrs=job_attrs,
                    shard_string=shard_string,
                    name_sort_for_merge=True,
                )
                # Name-sorting with samtools; deduplication runs after the merge, since dupblaster
                # must see the full RG-ordered read stream to catch cross-shard duplicates.
                if shard_fifo_epilogue:
                    cmd += f'\n{shard_fifo_epilogue}'
                j.command(cmd)
                sorted_bams.append(j.sorted_bam)
                sharded_align_jobs.append(j)

        merge_j = batch_instance.new_bash_job('Merge BAMs', (job_attrs or {}) | {'tool': 'samtools_merge'})
        merge_j.image(config.config_retrieve(['workflow', 'driver_image']))

        # for FASTQ or BAM inputs, requesting more disk (400G). Example when
        # default is not enough: https://batch.hail.populationgenomics.org.au/batches/73892/jobs/56
        nthreads = config.config_retrieve(['workflow', 'merge_threads'], 16)
        storage_gb = storage_for_align_job(alignment_input=alignment_input)
        merge_j.storage(f'{storage_gb}GiB')
        merge_j.cpu(nthreads)
        merge_j.memory('highmem')

        # Merge the name-sorted shards in name order, deduplicate the combined RG-ordered
        # stream with dupblaster, then coordinate-sort the result for downstream tools.
        align_cmd = f"""\
        {PIPEFAIL}
        samtools merge -n -@{min(nthreads, 6) - 1} - {' '.join(map(str, sorted_bams))} \
        """.strip()
        jobs.extend(sharded_align_jobs)
        jobs.append(merge_j)
        merge_or_align_j = merge_j

    jobs.append(
        finalise_alignment(
            align_cmd=align_cmd,
            job=merge_or_align_j,
            output_path=output_path,
            out_markdup_metrics_path=markdup_metrics_path,
            fifo_epilogue=fifo_epilogue,
        )
    )

    return jobs


def storage_for_align_job(alignment_input: filetypes.AlignmentInput) -> int:
    """
    Get storage for an alignment job, gb
    This function pulls the value from the config, otherwise utilises the full capacity of a standard
    machine. Increases storage if genomes are being handled, and again for unindexed/unsorted CRAM/BAMs.
    """

    if storage_gb := config.config_retrieve(['workflow', 'align_storage'], None):
        return storage_gb

    # Taking a full instance without attached by default:
    storage_gb = 100

    if config.config_retrieve(['workflow', 'sequencing_type']) == 'genome':
        storage_gb = 400

        # More disk is needed for FASTQ or BAM inputs than for realignment from CRAM
        if isinstance(alignment_input, filetypes.FastqPair | filetypes.FastqPairs | filetypes.BamPath):
            storage_gb = 400

        # For unindexed/unsorted CRAM or BAM inputs, extra storage is needed for tmp
        if isinstance(alignment_input, filetypes.CramPath | filetypes.BamPath) and not alignment_input.index_path:
            storage_gb += 150

    return storage_gb


def _align_one(
    alignment_input: filetypes.FastqPair | filetypes.FastqOraPair | filetypes.CramPath | filetypes.BamPath,
    sequencing_group_name: str,
    job_attrs: dict,
    shard_string: str | None = None,
    name_sort_for_merge: bool = False,
) -> tuple[BashJob, str, str]:
    """
    Creates a job that (re)aligns reads to hg38. Returns the job object and a command
    separately, and doesn't add the command to the Job object, so stream-sorting
    and/or deduplication can be appended to the command later.

    Note: When this function is called within the align function, DRAGMAP is used as the default tool.
    """

    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_bash_job(
        name=f'Align {shard_string} ' if shard_string else f'Align {alignment_input}',
        attributes=job_attrs | {'tool': 'dragmap'},
    )

    # allow alignment jobs to be non-spot
    job.spot(config.config_retrieve(['workflow', 'align_spot'], True))

    # stop messing about with the indirection of the resource profiles - set explicitly
    storage_gb = storage_for_align_job(alignment_input)
    nthreads = config.config_retrieve(['workflow', 'align_threads'], 16)
    # this uses the DRAGMAP base image, but with the ORAD decompression software added in (see config)
    job.storage(f'{storage_gb}GiB').cpu(nthreads).memory(
        config.config_retrieve(['workflow', 'align_memory'], 'highmem')
    ).image(config.config_retrieve(['images', 'dragmap']))

    sort_index_input_cmd = ''

    # 2022-07-22 mfranklin:
    #   Replace process substitution with named-pipes (FIFO)
    #   This is named-pipe name -> command to populate it
    fifo_commands: dict[str, str] = {}

    if isinstance(alignment_input, filetypes.CramPath | filetypes.BamPath):
        use_bazam = True
        shard_param = f' -s {shard_string.replace("/", ",")}' if shard_string else ''

        bazam_ref_cmd = ''
        samtools_ref_cmd = ''
        if isinstance(alignment_input, filetypes.CramPath):
            reference_inp = batch_instance.read_input_group(
                base=str(alignment_input.reference_assembly),
                fai=str(alignment_input.reference_assembly) + '.fai',
            ).base
            bazam_ref_cmd = f'-Dsamjdk.reference_fasta={reference_inp}'
            samtools_ref_cmd = f'--reference {reference_inp}'

        group = alignment_input.resource_group(batch_instance)

        if not alignment_input.index_path:
            sort_index_input_cmd = dedent(
                f"""
            mkdir -p $BATCH_TMPDIR/sorted
            mkdir -p $BATCH_TMPDIR/sort_tmp
            samtools sort {samtools_ref_cmd} \
            {group[alignment_input.ext]} \
            -@{nthreads - 1} \
            -T $BATCH_TMPDIR/sort_tmp \
            > $BATCH_TMPDIR/sorted.{alignment_input.ext}

            mv $BATCH_TMPDIR/sorted.{alignment_input.ext} {group[alignment_input.ext]}
            rm -rf $BATCH_TMPDIR/sort_tmp

            # bazam requires an index at foo.crai (not foo.cram.crai) so we must set the
            # path explicitly. We can not access the localized cram file's basename via the
            # Input ResourceGroup so using shell magic to strip the trailing "m" and add an
            # "i" to the alignment_path. This should work for both .cram and .bam files.
            alignment_path="{group[alignment_input.ext]}"
            samtools index -@{nthreads - 1} $alignment_path  ${{alignment_path%m}}i
            """,
            )

        bazam_cmd = dedent(
            f"""\
        bazam -Xmx16g {bazam_ref_cmd} \
        -n{min(nthreads, 6)} -bam {group[alignment_input.ext]}{shard_param} \
        """,
        )
        prepare_fastq_cmd = ''
        r1_param = 'r1'
        r2_param = ''
        fifo_commands[r1_param] = bazam_cmd

    # new functionality,
    elif isinstance(alignment_input, filetypes.FastqOraPair):
        fqo_resource_group = alignment_input.resource_group(batch_instance)
        use_bazam = False
        r1_param = '$BATCH_TMPDIR/R1.fq.gz'
        r2_param = '$BATCH_TMPDIR/R2.fq.gz'
        # Need file names to end with ".gz" for BWA or DRAGMAP to parse correctly:
        prepare_fastq_cmd = dedent(
            f"""
            {PIPEFAIL}
            tar -xf {fqo_resource_group.reference} -C $BATCH_TMPDIR
            orad -c --ora-reference $BATCH_TMPDIR/oradata_homo_sapiens {fqo_resource_group.r1} > {r1_param}
            rm {fqo_resource_group.r1}
            orad -c --ora-reference $BATCH_TMPDIR/oradata_homo_sapiens {fqo_resource_group.r2} > {r2_param}
            rm {fqo_resource_group.r2}
        """,
        )

    else:  # only for BAMs that are missing index
        fastq_pair = alignment_input.as_resources(batch_instance)
        use_bazam = False
        r1_param = '$BATCH_TMPDIR/R1.fq.gz'
        r2_param = '$BATCH_TMPDIR/R2.fq.gz'
        # Need file names to end with ".gz" for BWA or DRAGMAP to parse correctly:
        prepare_fastq_cmd = dedent(
            f"""\
        {PIPEFAIL}
        mv {fastq_pair.r1} {r1_param}
        mv {fastq_pair.r2} {r2_param}
        """,
        )

    dragmap_index = batch_instance.read_input_group(
        **{
            k.replace('.', '_'): os.path.join(config.config_retrieve(['references', 'dragmap_prefix']), k)
            for k in DRAGMAP_INDEX_FILES
        },
    )
    input_params = f'--interleaved=1 -b {r1_param}' if use_bazam else f'-1 {r1_param} -2 {r2_param}'
    # TODO: consider reverting to use of all threads if node capacity
    # issue is resolved: https://hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/Job.20becomes.20unresponsive
    # add a watch_disk loop to monitor disk space, and terminate the job before entering a zombie state.
    # default value is 2GB in KB
    storage_buffer = config.config_retrieve(['workflow', 'align_buffer_kb'], 2097152)
    cmd = f"""\
    watch_disk() {{
    local min_free_kb={storage_buffer}
    while true; do
      local avail
      avail=$(df --output=avail "$BATCH_TMPDIR" | tail -1)
      if (( avail < min_free_kb )); then
        echo "FATAL: $BATCH_TMPDIR has ${{avail}}KB free (< ${{min_free_kb}}KB) — aborting before disk fills" >&2
        kill -TERM -$$ 2>/dev/null
        exit 1
      fi
      sleep 15
    done
    }}
    watch_disk &
    WATCHDOG_PID=$!
    trap 'kill "$WATCHDOG_PID" 2>/dev/null' EXIT

    {prepare_fastq_cmd}
    dragen-os -r {dragmap_index} {input_params} \\
        --RGID {sequencing_group_name} --RGSM {sequencing_group_name} \\
        --num-threads {nthreads - 1}
    """

    # name-sort each shard so `samtools merge -n` produces an RG-ordered stream for dupblaster
    cmd = dedent(cmd).strip()
    if name_sort_for_merge:
        cmd += ' ' + name_sort_cmd(nthreads) + f' -o {job.sorted_bam}'

    # Wait-check appended after the core command by the caller, once any downstream pipe
    # stages (e.g. dedup/sort) have been attached - it must not sit between the aligner's
    # stdout and those stages, or the pipe ends up reading the echo output instead.
    fifo_epilogue = ''
    if fifo_commands:
        fifo_pre = [
            dedent(
                f"""
            mkfifo {fname}
            {sub_cmd} > {fname} &
            pid_{fname}=$!
        """,
            ).strip()
            for fname, sub_cmd in fifo_commands.items()
        ]

        _fifo_waits = ' && '.join(f'wait $pid_{fname}' for fname in fifo_commands)
        fifo_epilogue = dedent(
            f"""
            if {_fifo_waits}
            then
                echo -e "Background processes finished successfully"
            else
                # Background processes failed
                trap 'error1' ERR
            fi
        """,
        ).strip()

        # Now prepare command
        cmd = '\n'.join([PIPEFAIL, sort_index_input_cmd, *fifo_pre, cmd])
    return job, cmd, fifo_epilogue


def name_sort_cmd(nthreads: int) -> str:
    """
    Create command that name-sorts a SAM/BAM stream, grouping reads by name (RG order)
    so that dupblaster can deduplicate in-stream after the shards are merged.
    """
    return f'| samtools sort -n -@{min(nthreads, 6) - 1} -T $BATCH_TMPDIR/sam-sort-tmp -Obam '


def dedup_sort_cmd(nthreads: int, stats_path: str) -> str:
    """
    Create command that deduplicates a name-sorted (RG-ordered) stream with dupblaster,
    writing duplication stats to `stats_path`, then coordinate-sorts the result.
    """
    cmd = f'| dupblaster --stats {stats_path} '
    cmd += f'| samtools sort -@{min(nthreads, 6) - 1} -T $BATCH_TMPDIR/samtools-dd-tmp -Obam '
    return cmd


def finalise_alignment(
    align_cmd: str,
    job: BashJob,
    output_path: Path,
    out_markdup_metrics_path: str,
    fifo_epilogue: str = '',
) -> BashJob:
    """
    Deduplication (dupblaster) and CRAM conversion (samtools) both run in-stream within the
    alignment/merge job - the DRAGMAP and samtools images both provide samtools, so no extra
    job or image is needed. This appends the CRAM conversion to the pipeline and persists the
    duplication stats and indexed CRAM.
    """

    batch_instance = hail_batch.get_batch()

    nthreads = config.config_retrieve(['workflow', 'merge_threads'], 8)

    job.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        },
    )
    fasta_reference = hail_batch.fasta_res_group(batch_instance)

    align_cmd = align_cmd.strip()

    # dedup the raw RG-ordered dragen-os stream, then coordinate-sort
    align_cmd += f' {dedup_sort_cmd(nthreads, job.markdup_metrics)}'

    # convert the coordinate-sorted, duplicate-marked stream to an indexed CRAM in-stream
    align_cmd += (
        f' | samtools view --write-index -@{min(nthreads, 6) - 1} '
        f'-T {fasta_reference.base} -O cram,version=3.0 -o {job.output_cram.cram} -'
    )
    if fifo_epilogue:
        # bazam background-writer wait-check: must come after the full pipeline, not between
        # the aligner and the dedup/sort/view stages consuming its stdout
        align_cmd += f'\n{fifo_epilogue}'

    job.command(align_cmd)

    # persist the duplication stats produced by dupblaster and the indexed CRAM
    batch_instance.write_output(job.markdup_metrics, out_markdup_metrics_path)
    batch_instance.write_output(job.output_cram, str(output_path.with_suffix('')))

    return job
