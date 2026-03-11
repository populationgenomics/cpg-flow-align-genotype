import argparse
from cpg_utils import hail_batch
from hailtop.batch.job import Job

def run_fixmate(batch: hail_batch.Batch, bam_path: str, out_bam_path: str) -> Job:
    job = batch.new_job('Fix Primary Alignments with Validation')

    job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/cpg-flow-align-genotype:0.4.5-2')
    job.memory('16Gi')
    job.storage('250Gi')

    cpu_count = 8
    job.cpu(cpu_count)

    bam_localised = batch.read_input_group(bam=bam_path, bai=f'{bam_path}.bai')

    job.declare_resource_group(
        output_bam={'bam': '{root}.bam', 'bam.bai': '{root}.bam.bai'}
    )

    job.command(f'''
        set -ex
        
        # 1. Sort by NAME to group duplicates together
        samtools sort -n -@ {cpu_count} -o name_sorted.bam {bam_localised.bam} -T $BATCH_TMPDIR

        # 2. Process with Python
        cat > fix_primary.py << EOF
import pysam

last_name = None
r1_primary_found = False
r2_primary_found = False
modified_count = 0
total_count = 0

with pysam.AlignmentFile("name_sorted.bam", "rb") as inf:
    with pysam.AlignmentFile("fixed_name_sorted.bam", "wb", template=inf) as outf:
        for read in inf:
            total_count += 1
            
            # If we hit a new read name, reset our primary trackers
            if read.query_name != last_name:
                last_name = read.query_name
                r1_primary_found = False
                r2_primary_found = False

            # We only care about reads not already marked secondary/supplementary
            if not (read.is_secondary or read.is_supplementary):
                if read.is_read1:
                    if r1_primary_found:
                        read.is_supplementary = True
                        modified_count += 1
                    else:
                        r1_primary_found = True
                else:
                    if r2_primary_found:
                        read.is_supplementary = True
                        modified_count += 1
                    else:
                        r2_primary_found = True

            outf.write(read)

        print(f"--- Processing Summary ---")
        print(f"Total reads processed: {{total_count}}")
        print(f"Reads converted to supplementary: {{modified_count}}")
EOF
        
        # Run the Python script using the Python interpreter
        /usr/bin/env python fix_primary.py
        
        # 3. Restore Coordinate Sort order and Index
        samtools sort -@ {cpu_count} -o {job.output_bam.bam} fixed_name_sorted.bam -T $BATCH_TMPDIR
        samtools index {job.output_bam.bam}
    ''')

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

