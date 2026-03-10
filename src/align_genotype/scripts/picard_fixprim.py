import argparse
from cpg_utils import config, hail_batch
from hailtop.batch.job import Job

def run_fixmate(batch: hail_batch.Batch, bam_path: str, out_bam_path: str) -> Job:
    """Run Picard FixMateInformation on BAM file to clean mate information."""

    job = batch.new_job(
        'Picard FixMateInformation',
    )

    # Use the custom image that has pysam installed
    job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/cpg-flow-align-genotype:0.4.5-2')
    job.memory('16Gi')
    job.storage('250Gi')
    job.cpu(4)

    # read in the BAM and index
    bam_localised = batch.read_input_group(
        bam=bam_path,
        bai=f'{bam_path}.bai',
    ).bam

    # Declare output resource group
    job.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai',
        },
    )

    # Use a set to track read names we've already "assigned" as primary
    job.command(f'''
        set -ex
        
        cat > fix_primary.py << 'EOF'
import pysam

seen_primary_r1 = set()
seen_primary_r2 = set()

with pysam.AlignmentFile("{bam_localised}", "rb") as inf:
    with pysam.AlignmentFile("{job.output_bam.bam}", "wb", template=inf) as outf:
        for read in inf:
            # Skip reads that are already secondary or supplementary
            if read.is_secondary or read.is_supplementary:
                outf.write(read)
                continue

            # Check if this is Read 1 or Read 2
            current_set = seen_primary_r1 if read.is_read1 else seen_primary_r2

            if read.query_name in current_set:
                # We have already seen a primary for this name/pair-end.
                # Convert this one to SUPPLEMENTARY (0x800)
                read_type = "R1" if read.is_read1 else "R2"
                print(f"Changed {{read.query_name}} ({{read_type}}) from primary to supplementary")
                read.is_supplementary = True
            else:
                # First time seeing this name, keep it as Primary
                current_set.add(read.query_name)

            outf.write(read)
EOF
        
        # Run the Python script using the Python interpreter
        /usr/bin/env python fix_primary.py
        
        # Index the output BAM
        samtools index {job.output_bam.bam}
    '''
    )

    batch.write_output(job.output_bam, out_bam_path.removesuffix('.bam'))
    return job

def main(args):
    batch = hail_batch.get_batch()

    run_fixmate(batch, bam_path=args.bam_path, out_bam_path=args.out_bam_path)

    batch.run(wait=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_path', type=str, required=True, help='Path to input BAM file')
    parser.add_argument('--out_bam_path', type=str, required=True, help='Path to output BAM file')
    args = parser.parse_args()
    main(args)

