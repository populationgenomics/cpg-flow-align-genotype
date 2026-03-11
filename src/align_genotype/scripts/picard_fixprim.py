import argparse
from cpg_utils import hail_batch
from hailtop.batch.job import Job

def job_name_sort(batch: hail_batch.Batch, bam_path: str) -> Job:
    """Job 1: Sort BAM by name to group duplicate alignments together."""
    job = batch.new_job('Name Sort BAM')

    job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/cpg-flow-align-genotype:0.4.5-2')
    job.memory('16Gi')
    job.storage('500Gi')

    cpu_count = 8
    job.cpu(cpu_count)

    bam_localised = batch.read_input_group(bam=bam_path, bai=f'{bam_path}.bai')

    job.declare_resource_group(
        name_sorted={'bam': '{root}.bam'}
    )

    job.command(f'''
        set -ex
        samtools sort -n -@ {cpu_count} -o {job.name_sorted.bam} {bam_localised.bam} -T $BATCH_TMPDIR
    ''')

    return job

def job_fix_primary_alignments(batch: hail_batch.Batch, name_sorted_bam) -> Job:
    """Job 2: Process name-sorted BAM to fix duplicate primary alignments."""
    job = batch.new_job('Fix Primary Alignments')

    job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/cpg-flow-align-genotype:0.4.5-2')
    job.memory('16Gi')
    job.storage('500Gi')
    job.cpu(4)

    job.declare_resource_group(
        fixed={'bam': '{root}.bam'}
    )

    job.command(f'''
        set -ex
        
        cat > fix_primary.py << EOF
import pysam
import random

def get_alignment_score(read):
    """Get alignment score from AS tag, or use MAPQ if AS not present."""
    try:
        return read.get_tag('AS')
    except KeyError:
        return None

def get_sort_key(read):
    """Return tuple for sorting: (MAPQ, AS, random). Higher is better."""
    as_score = get_alignment_score(read)
    # Use -1 for missing AS so it sorts lower than any real score
    return (read.mapping_quality, as_score if as_score is not None else -1, random.random())

last_name = None
current_r1_reads = []
current_r2_reads = []
modified_count = 0
total_count = 0

def process_read_group(reads, read_type):
    """Process a group of reads with the same name and type (R1 or R2)."""
    global modified_count
    
    if len(reads) <= 1:
        return reads
    
    # Filter to only primary alignments (not already secondary/supplementary)
    primary_reads = [r for r in reads if not (r.is_secondary or r.is_supplementary)]
    
    if len(primary_reads) <= 1:
        return reads
    
    # Sort by MAPQ (desc), AS (desc), random
    # The best one will be at the end after sorting
    sorted_reads = sorted(primary_reads, key=get_sort_key)
    best_read = sorted_reads[-1]
    
    # Mark all others as supplementary
    result_reads = []
    for read in reads:
        if read is best_read or read.is_secondary or read.is_supplementary:
            result_reads.append(read)
        else:
            # This read needs to be marked as supplementary
            mapq = read.mapping_quality
            as_score = get_alignment_score(read)
            as_str = "AS=" + str(as_score) if as_score is not None else "AS=None"
            chrom = read.reference_name if read.reference_name else "unmapped"
            pos = str(read.reference_start) if read.reference_start is not None else "N/A"
            
            print("BEFORE: " + read.query_name + " (" + read_type + ") - " + chrom + ":" + pos + " MAPQ=" + str(mapq) + " " + as_str + " Flags=" + str(read.flag))
            
            read.is_supplementary = True
            modified_count += 1
            
            print("AFTER:  " + read.query_name + " (" + read_type + ") - " + chrom + ":" + pos + " MAPQ=" + str(mapq) + " " + as_str + " Flags=" + str(read.flag) + " [NOW SUPPLEMENTARY]")
            print()
            
            result_reads.append(read)
    
    return result_reads

with pysam.AlignmentFile("{name_sorted_bam}", "rb") as inf:
    with pysam.AlignmentFile("{job.fixed.bam}", "wb", template=inf) as outf:
        for read in inf:
            total_count += 1
            
            # If we hit a new read name, process the accumulated reads
            if read.query_name != last_name:
                if last_name is not None:
                    # Process previous read group
                    for r in process_read_group(current_r1_reads, "R1"):
                        outf.write(r)
                    for r in process_read_group(current_r2_reads, "R2"):
                        outf.write(r)
                
                # Reset for new read name
                last_name = read.query_name
                current_r1_reads = []
                current_r2_reads = []
            
            # Accumulate reads by type
            if read.is_read1:
                current_r1_reads.append(read)
            else:
                current_r2_reads.append(read)
        
        # Don't forget the last group
        if last_name is not None:
            for r in process_read_group(current_r1_reads, "R1"):
                outf.write(r)
            for r in process_read_group(current_r2_reads, "R2"):
                outf.write(r)

print("--- Processing Summary ---")
print("Total reads processed: " + str(total_count))
print("Reads converted to supplementary: " + str(modified_count))
EOF
        
        /usr/bin/env python fix_primary.py
    ''')

    return job

def job_coordinate_sort_and_index(batch: hail_batch.Batch, fixed_bam, out_bam_path: str) -> Job:
    """Job 3: Sort BAM by coordinate and create index."""
    job = batch.new_job('Coordinate Sort and Index')

    job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/cpg-flow-align-genotype:0.4.5-2')
    job.memory('16Gi')
    job.storage('500Gi')

    cpu_count = 8
    job.cpu(cpu_count)

    job.declare_resource_group(
        output_bam={'bam': '{root}.bam', 'bam.bai': '{root}.bam.bai'}
    )

    job.command(f'''
        set -ex
        samtools sort -@ {cpu_count} -o {job.output_bam.bam} {fixed_bam} -T $BATCH_TMPDIR
        samtools index {job.output_bam.bam}
    ''')

    batch.write_output(job.output_bam, out_bam_path.removesuffix('.bam'))
    return job

def main(args):
    batch = hail_batch.get_batch()

    # Job 1: Name sort
    job1 = job_name_sort(batch, bam_path=args.bam_path)

    # Job 2: Fix primary alignments (depends on job1)
    job2 = job_fix_primary_alignments(batch, name_sorted_bam=job1.name_sorted.bam)

    # Job 3: Coordinate sort and index (depends on job2)
    job3 = job_coordinate_sort_and_index(batch, fixed_bam=job2.fixed.bam, out_bam_path=args.out_bam_path)

    batch.run(wait=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_path', type=str, required=True, help='Path to input BAM file')
    parser.add_argument('--out_bam_path', type=str, required=True, help='Path to output BAM file')
    args = parser.parse_args()
    main(args)

