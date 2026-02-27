"""
Intention here is to run a variety of tools for duplicate marking, and evaluate their runtime performance
Eventually this will lead to a proper evaluation, comparing the resulting alignments

Making some big storage/cpu/memory assumptions, and can scale all of that back down once we observe performance
"""

from argparse import ArgumentParser
from hailtop.batch import Batch
from hailtop.batch.job import BashJob

from cpg_utils import config, hail_batch
from cpg_flow import utils as flow_utils


# I'm embedding this directly to make target images super obvious
IMAGES: dict[str, str] = {
    'rust_dupmark': 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/rust_dupmarker:0.1.0-1',
    'samblaster': 'australia-southeast1-docker.pkg.dev/cpg-common/images/samblaster:0.1.26-1',
    'sambamba': 'australia-southeast1-docker.pkg.dev/cpg-common/images/sambamba:1.0.1-1',
    'streammd': 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/streammd:4.3.0-1',
}


def make_a_job(batch: Batch, tool: str) -> BashJob:
    """ronseal."""
    new_job = batch.new_bash_job(name=f'Tool: {tool}')
    new_job.image(IMAGES[tool])
    new_job.cpu(8)
    new_job.memory('highmem')
    new_job.storage('360GiB')
    return new_job


def create_sambamba_job(batch: Batch, bamfile: str, outfile: str) -> None:
    """
    https://lomereiter.github.io/sambamba/docs/sambamba-markdup.html
    sambamba markdup OPTIONS <input.bam> <output.bam>

    only works BAM -> BAM, needs further CRAM compression step
    """

    if flow_utils.exists(outfile):
        return

    sb_job = make_a_job(batch, 'sambamba')

    sb_job.command(f"""
    sambamba markdup \\
        -t 4 \\
        --tmpdir=$BATCH_TMPDIR \\
        {bamfile} {sb_job.out}
    """)
    batch.write_output(sb_job.out, outfile)


def create_samblaster_job(batch: Batch, bamfile: str, reference: str):
    """ronseal."""
    sb_job = make_a_job(batch, 'samblaster')


def create_streammd_job(batch: Batch, bamfile: str, outfile: str):
    """
    e.g. bwa mem ref.fa r1.fq r2.fq|streammd
    """

    if flow_utils.exists(outfile):
        return

    streammd_job = make_a_job(batch, 'streammd')
    streammd_job.command(f"""
    samtools view \\
        -h \\
        --output-fmt SAM \\
        {bamfile} | \\
    streammd \
        --output {streammd_job.out}
    """)

    batch.write_output(streammd_job.out, outfile)


def create_dupmark_job(batch: Batch, bamfile: str, reference: str, outfile: str):
    """
    A high-performance SAM/BAM/CRAM duplicate marker using a Bloom filter.
    Reads from standard input and writes CRAM to standard output

    Usage: dupmark [OPTIONS]

    Options:
        -e, --expected-items <EXPECTED_ITEMS>
        Expected number of unique reads to process. This sizes the bloom filter [default: 100000000]
        -f, --false-pos-rate <FALSE_POS_RATE>
        Acceptable false positive rate for the bloom filter [default: 0.0001]
        -r, --reference <REFERENCE>
        Optional reference FASTA for CRAM compression
    """

    if flow_utils.exists(outfile):
        return

    dm_job = make_a_job(batch, 'rust_dupmark')
    dm_job.command(f"""
    cat {bamfile} | dupmark -r {reference} > {dm_job.out}
    """)
    batch.write_output(dm_job.out, outfile)


def main(bamfile: str, outdir: str) -> None:
    batch_instance = hail_batch.get_batch(name='Various Duplicate Markers test.')

    # set with workflow.ref_fasta
    ref_fa = hail_batch.fasta_res_group(batch_instance).base

    input_bam = batch_instance.read_input(bamfile)

    create_sambamba_job(batch_instance, input_bam, outfile=f'{outdir}/sambamba/result.bam')

    create_streammd_job(batch_instance, input_bam, outfile=f'{outdir}/streammd/result.bam')

    create_dupmark_job(batch_instance, input_bam, outfile=f'{outdir}/rust_dupmark/result.cram', reference=ref_fa)

    batch_instance.run(wait=False)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='Sorted BAM file', required=True)
    parser.add_argument('-o', help='output root, all outputs derived from this path', required=True)
    args = parser.parse_args()
    main(bamfile=args.i, outdir=args.o)
