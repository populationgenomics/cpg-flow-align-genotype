"""
NOT DONE YET

Secondary part of the process - now that there are some duplicate-marked files (mix of CRAM and BAM),
run them through the variant calling process we already have, to create a gVCF per dupmark tool
"""

from cpg_utils import hail_batch, to_path

from align_genotype.jobs.genotype import genotype


CRAMS = {
    'rust_dupmark': 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/cram30/rust_dupmark/result.cram',
    'sambamba': 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/cram30/sambamba/result.cram',
    'samblaster': 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/samblaster/result.cram',
    # 'streammd': 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/streammd/result.cram.cram',  # not completed yet
}


def main() -> None:
    batch_instance = hail_batch.get_batch(name='Various Duplicate Markers test.')

    # get the intervals file from main workflow outputs (copied over)
    intervals_path = 'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/intervals'
    intervals_as_paths = list(to_path(intervals_path).glob('*.interval_list'))

    for marker, cramfile in CRAMS.items():
        _jobs = genotype(
            sequencing_group_name=marker,
            tmp_prefix=to_path(f'gs://cpg-ghfm-kidgen-test-tmp/duplicate_marker_test'),
            cram_path=cramfile,
            output_path=to_path(f'gs://cpg-ghfm-kidgen-test/duplicate_marker_test/haplotypecaller/{marker}.g.vcf.gz'),
            intervals=intervals_as_paths,
            job_attrs={'tool': marker}
        )

    batch_instance.run(wait=False)


if __name__ == '__main__':
    main()
