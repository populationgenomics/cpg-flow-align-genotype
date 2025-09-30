"""
FASTQ/BAM/CRAM -> CRAM: create Hail Batch jobs for (re-)alignment.
"""

import logging
import os.path
from textwrap import dedent
from typing import cast

from cpg_flow import filetypes, resources, targets, utils
from cpg_utils import Path, config, hail_batch, to_path
from hailtop.batch.job import BashJob

from align_genotype.jobs import picard

DRAGMAP_INDEX_FILES = ['hash_table.cfg.bin', 'hash_table.cmp', 'reference.bin']


def align(
    sequencing_group: targets.SequencingGroup,
    job_attrs: dict,
    output_path: Path,
    sorted_bam_path: str,
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

    - if the sorted bam can be reused, skip the alignment job(s) and go straight to markdup.

    - the markdup tool is picard, deduplication is submitted in a separate job.

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

    if utils.can_reuse(sorted_bam_path):
        logging.info(f'{sequencing_group.id} :: Re-using sorted BAM: {sorted_bam_path}')
        # It's necessary to create this merge_or_align_j object to pass it to finalise_alignment,
        # and to declare the sorted_bam_path as a resource, so it can be written to the checkpoint.
        job_attrs = (job_attrs or {}) | {'label': 'Reusing sorted bam', 'tool': 'Reusing sorted bam'}
        merge_or_align_j = batch_instance.new_bash_job('Reusing sorted bam', job_attrs or {})
        merge_or_align_j.image(config.config_retrieve(['workflow', 'driver_image']))
        merge_or_align_j.sorted_bam = batch_instance.read_input(sorted_bam_path)
        jobs.append(merge_or_align_j)
        # The align_cmd and other parameters are not used but are necessary to pass to finalise_alignment.
        align_cmd = ''
        stdout_is_sorted = True

    elif not sharded:  # Just running one alignment job
        if isinstance(alignment_input, filetypes.FastqPairs):
            alignment_input = alignment_input[0]
        align_j, align_cmd = _align_one(
            alignment_input=alignment_input,
            sequencing_group_name=sequencing_group.id,
            job_attrs=job_attrs,
            should_sort=False,
        )
        stdout_is_sorted = False
        jobs.append(align_j)
        merge_or_align_j = align_j

    # Aligning in parallel and merging afterwards
    else:
        if sharded_fq:  # Aligning each lane separately, merging after
            # running alignment for each fastq pair in parallel
            fastq_pairs = cast('filetypes.FastqPairs', alignment_input)
            for pair in fastq_pairs:
                # dragmap command, without sorting and deduplication:
                j, cmd = _align_one(
                    alignment_input=pair,
                    sequencing_group_name=sequencing_group.id,
                    job_attrs=job_attrs,
                    should_sort=True,
                )
                j.command(hail_batch.command(cmd, monitor_space=True))
                sorted_bams.append(j.sorted_bam)
                sharded_align_jobs.append(j)

        elif sharded_bazam:  # Using BAZAM to shard CRAM
            for shard_number in range(realignment_shards):
                shard_string = f'{shard_number + 1}/{realignment_shards}'
                j, cmd = _align_one(
                    alignment_input=alignment_input,
                    sequencing_group_name=sequencing_group.id,
                    job_attrs=job_attrs,
                    shard_string=shard_string,
                    should_sort=True,
                )
                # Sorting with samtools, but not adding deduplication yet, because we need to merge first.
                j.command(hail_batch.command(cmd, monitor_space=True))
                sorted_bams.append(j.sorted_bam)
                sharded_align_jobs.append(j)

        merge_j = batch_instance.new_bash_job('Merge BAMs', (job_attrs or {}) | {'tool': 'samtools_merge'})
        merge_j.image(config.config_retrieve(['images', 'samtools']))

        nthreads = resources.STANDARD.set_resources(
            j=merge_j,
            nthreads=config.config_retrieve(['workflow', 'align_threads'], resources.STANDARD.max_threads()),
            # for FASTQ or BAM inputs, requesting more disk (400G). Example when
            # default is not enough: https://batch.hail.populationgenomics.org.au/batches/73892/jobs/56
            storage_gb=storage_for_align_job(
                alignment_input=alignment_input,
            ),
        ).get_nthreads()

        align_cmd = f"""\
        samtools merge -@{nthreads - 1} - {' '.join(map(str, sorted_bams))}
        """.strip()
        jobs.extend(sharded_align_jobs)
        jobs.append(merge_j)
        merge_or_align_j = merge_j
        stdout_is_sorted = True

    jobs.append(
        finalise_alignment(
            align_cmd=align_cmd,
            stdout_is_sorted=stdout_is_sorted,
            job=merge_or_align_j,
            job_attrs=job_attrs,
            output_path=output_path,
            sorted_bam_path=sorted_bam_path,
            out_markdup_metrics_path=markdup_metrics_path,
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
    storage_gb = resources.STANDARD.calc_instance_disk_gb()
    if config.config_retrieve(['workflow', 'sequencing_type']) == 'genome':
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
    should_sort: bool = False,
) -> tuple[BashJob, str]:
    """
    Creates a job that (re)aligns reads to hg38. Returns the job object and a command
    separately, and doesn't add the command to the Job object, so stream-sorting
    and/or deduplication can be appended to the command later.

    Note: When this function is called within the align function, DRAGMAP is used as the default tool.

    Args:
        batch_instance ():
        alignment_input ():
        sequencing_group_name ():
        job_attrs ():
        shard_string (str | None): if populated, a String in the form 'X/Y', this is string formatted as a CLI argument
        should_sort ():
    """

    batch_instance = hail_batch.get_batch()

    requested_nthreads = config.config_retrieve(['workflow', 'align_threads'], resources.STANDARD.max_threads())

    job_name = f'Align {shard_string} ' if shard_string else f'Align {alignment_input}'

    job = batch_instance.new_bash_job(name=job_name, attributes=job_attrs | {'tool': 'dragmap'})

    vm_type = resources.HIGHMEM if config.config_retrieve(['workflow', 'align_use_highmem']) else resources.STANDARD

    # confused by this error message
    nthreads = vm_type.set_resources(
        j=job,
        nthreads=config.config_retrieve(['workflow', 'align_threads'], resources.STANDARD.max_threads()),
        storage_gb=storage_for_align_job(alignment_input=alignment_input),
    ).get_nthreads()

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
        fastq_pair = alignment_input.as_resources(batch_instance)
        use_bazam = False
        r1_param = '$BATCH_TMPDIR/R1.fq.gz'
        r2_param = '$BATCH_TMPDIR/R2.fq.gz'
        # Need file names to end with ".gz" for BWA or DRAGMAP to parse correctly:
        prepare_fastq_cmd = dedent(
            f"""
            tar -xf {fastq_pair.reference} -C $BATCH_TMPDIR
            ls $BATCH_TMPDIR
            ls $BATCH_TMPDIR/**
            orad -c --ora-reference $BATCH_TMPDIR/oradata_homo_sapiens {fastq_pair.r1} > {r1_param}
            rm fastq_pair.r1
            orad -c --ora-reference $BATCH_TMPDIR/oradata_homo_sapiens {fastq_pair.r2} {r2_param}
            rm fastq_pair.r2
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
        mv {fastq_pair.r1} {r1_param}
        mv {fastq_pair.r2} {r2_param}
        """,
        )

    # this uses the DRAGMAP base image, but with the ORAD decompression software added in (see config)
    job.image(config.config_retrieve(['images', 'dragmap']))

    dragmap_index = batch_instance.read_input_group(
        **{
            k.replace('.', '_'): os.path.join(config.config_retrieve(['references', 'dragmap_prefix']), k)
            for k in DRAGMAP_INDEX_FILES
        },
    )
    input_params = f'--interleaved=1 -b {r1_param}' if use_bazam else f'-1 {r1_param} -2 {r2_param}'
    # TODO: consider reverting to use of all threads if node capacity
    # issue is resolved: https://hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/Job.20becomes.20unresponsive
    cmd = f"""\
    {prepare_fastq_cmd}
    dragen-os -r {dragmap_index} {input_params} \\
        --RGID {sequencing_group_name} --RGSM {sequencing_group_name} \\
        --num-threads {nthreads - 1}
    """

    # prepare command for adding sort on the end
    cmd = dedent(cmd).strip()
    if should_sort:
        cmd += ' ' + sort_cmd(requested_nthreads) + f' -o {job.sorted_bam}'

    if fifo_commands:
        fifo_pre = [
            dedent(
                f"""
            mkfifo {fname}
            {cmd} > {fname} &
            pid_{fname}=$!
        """,
            ).strip()
            for fname, cmd in fifo_commands.items()
        ]

        _fifo_waits = ' && '.join(f'wait $pid_{fname}' for fname in fifo_commands)
        fifo_post = dedent(
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
        cmd = '\n'.join([sort_index_input_cmd, *fifo_pre, cmd, fifo_post])
    return job, cmd


def sort_cmd(requested_nthreads: int) -> str:
    """
    Create command that coordinate-sorts SAM file
    """
    nthreads = resources.STANDARD.request_resources(nthreads=requested_nthreads).get_nthreads()
    return dedent(
        f"""\
    | samtools sort -@{min(nthreads, 6) - 1} -T $BATCH_TMPDIR/samtools-sort-tmp -Obam
    """,
    ).strip()


def finalise_alignment(
    align_cmd: str,
    stdout_is_sorted: bool,
    job: BashJob,
    job_attrs: dict,
    output_path: Path,
    sorted_bam_path: str,
    out_markdup_metrics_path: str,
) -> BashJob:
    """
    For `MarkDupTool.PICARD`, creates a new job, as Picard can't read from stdin.
    """

    batch_instance = hail_batch.get_batch()

    nthreads = resources.STANDARD.request_resources(
        nthreads=config.config_retrieve(
            ['workflow', 'align_threads'],
            resources.STANDARD.max_threads(),
        )
    ).get_nthreads()

    align_cmd = align_cmd.strip()
    if not stdout_is_sorted:
        align_cmd += f' {sort_cmd(nthreads)}'
    align_cmd += f' > {job.sorted_bam}'

    if not utils.can_reuse(sorted_bam_path):
        job.command(hail_batch.command(align_cmd, monitor_space=True))

    if (
        sorted_bam_path
        and not to_path(sorted_bam_path).exists()
        and config.config_retrieve(['workflow', 'checkpoint_sorted_bam'], False)
    ):
        # Write the sorted BAM to the checkpoint if it doesn't already exist and the config is set
        logging.info(f'Will write sorted bam to checkpoint: {sorted_bam_path}')
        batch_instance.write_output(job.sorted_bam, sorted_bam_path)

    return picard.markdup(
        batch_instance=batch_instance,
        sorted_bam=job.sorted_bam,
        output_path=output_path,
        out_markdup_metrics_path=out_markdup_metrics_path,
        job_attrs=job_attrs,
    )
