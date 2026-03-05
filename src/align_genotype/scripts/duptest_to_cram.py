"""
Rinse the various output files from duplicate marking into definitely-3.0-CRAM files
"""

from argparse import ArgumentParser
from hailtop.batch import Batch
from hailtop.batch.job import BashJob

from cpg_utils import hail_batch, to_path


# I'm embedding this directly to make target images super obvious
IMAGES: dict[str, str] = {
    'rust_dupmark': 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/rust_dupmarker:0.1.0-1',
    'samblaster': 'australia-southeast1-docker.pkg.dev/cpg-common/images/samblaster:0.1.26-1',
    'sambamba': 'australia-southeast1-docker.pkg.dev/cpg-common/images/sambamba:1.0.1-1',
    'streammd': 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/streammd:4.3.0-1',
    'samtools': 'australia-southeast1-docker.pkg.dev/cpg-common/images/samtools:1.21-1',
}

OUTPUTS: dict[str, str] = {
    'rust_dupmark': 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/rust_dupmark/result.cram',
    'sambamba': 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/sambamba/result.bam',
    # 'samblaster': 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/samblaster/result.cram',  # already handled
    # 'streammd': 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/streammd/result.cram',  # already handled
}

needs_sorting = ['samblaster', 'streammd']


def make_a_job(batch: Batch, tool: str) -> BashJob:
    """ronseal."""
    new_job = batch.new_bash_job(name=f'Convert: {tool}')
    new_job.image(IMAGES['samtools'])
    new_job.cpu(2)
    new_job.memory('highmem')
    new_job.storage('100GiB')
    return new_job


def main(outdir: str) -> None:
    batch_instance = hail_batch.get_batch(name='Duplicate Marker outputs to CRAM.')

    # set with workflow.ref_fasta
    ref_fa = hail_batch.fasta_res_group(batch_instance).base

    for tool, file in OUTPUTS.items():
        outroot = f'{outdir}/{tool}/result'

        if to_path(f'{outroot}.cram').exists():
            continue

        input_file = batch_instance.read_input(file)

        job = make_a_job(batch_instance, tool)

        job.declare_resource_group(
            output={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            },
        )

        job.command(f"""
        samtools view --write-index -@2 \\
            -T {ref_fa} \\
            -O cram,version=3.0 \\
            -o {job.output.cram} \\
            {input_file}
        echo "samtools view finished successfully"
        """)

        batch_instance.write_output(job.output, outroot)

    batch_instance.run(wait=False)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-o', help='Output root, all outputs derived from this path', required=True)
    args = parser.parse_args()
    main(outdir=args.o)
