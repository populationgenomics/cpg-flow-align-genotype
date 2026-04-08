"""
Script to fix duplicate primary alignments in a name-sorted BAM file.

This script processes a name-sorted BAM file and converts duplicate primary
alignments to supplementary alignments, keeping only the best alignment
(by MAPQ, then AS score) as primary.

Usage:
    python fix_primary_alignments.py --input <name_sorted.bam> --output <fixed.bam>
"""

import argparse
import random

import pysam


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
    return read.mapping_quality, as_score if as_score is not None else -1, random.random()  # noqa: S311


def process_read_group(reads, read_type, modified_count):
    """Process a group of reads with the same name and type (R1 or R2)."""
    if len(reads) <= 1:
        return reads, modified_count

    # Filter to only primary alignments (not already secondary/supplementary)
    primary_reads = [r for r in reads if not (r.is_secondary or r.is_supplementary)]

    if len(primary_reads) <= 1:
        return reads, modified_count

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
            as_str = 'AS=' + str(as_score) if as_score is not None else 'AS=None'
            chrom = read.reference_name if read.reference_name else 'unmapped'
            pos = str(read.reference_start) if read.reference_start is not None else 'N/A'

            print(f'BEFORE: {read.query_name} ({read_type}) - {chrom}:{pos} MAPQ={mapq} {as_str} Flags={read.flag}')

            read.is_supplementary = True
            modified_count += 1

            print(
                f'AFTER:  {read.query_name} ({read_type}) - {chrom}:{pos} '
                f'MAPQ={mapq} {as_str} Flags={read.flag} [NOW SUPPLEMENTARY]'
            )
            print()

            result_reads.append(read)

    return result_reads, modified_count


def fix_primary_alignments(input_bam: str, output_bam: str):
    """
    Process a name-sorted BAM file to fix duplicate primary alignments.

    Args:
        input_bam: Path to name-sorted input BAM file
        output_bam: Path to output BAM file with fixed alignments
    """
    last_name = None
    current_r1_reads: list[pysam.AlignedSegment] = []
    current_r2_reads: list[pysam.AlignedSegment] = []
    modified_count = 0
    total_count = 0

    with pysam.AlignmentFile(input_bam, 'rb') as inf, pysam.AlignmentFile(output_bam, 'wb', template=inf) as outf:
        for read in inf:
            total_count += 1

            # If we hit a new read name, process the accumulated reads
            if read.query_name != last_name:
                if last_name is not None:
                    # Process previous read group
                    processed_r1, modified_count = process_read_group(current_r1_reads, 'R1', modified_count)
                    for r in processed_r1:
                        outf.write(r)
                    processed_r2, modified_count = process_read_group(current_r2_reads, 'R2', modified_count)
                    for r in processed_r2:
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
            processed_r1, modified_count = process_read_group(current_r1_reads, 'R1', modified_count)
            for r in processed_r1:
                outf.write(r)
            processed_r2, modified_count = process_read_group(current_r2_reads, 'R2', modified_count)
            for r in processed_r2:
                outf.write(r)

    print('--- Processing Summary ---')
    print(f'Total reads processed: {total_count}')
    print(f'Reads converted to supplementary: {modified_count}')


def main():
    parser = argparse.ArgumentParser(description='Fix duplicate primary alignments in a name-sorted BAM file')
    parser.add_argument('--input', type=str, required=True, help='Path to name-sorted input BAM file')
    parser.add_argument('--output', type=str, required=True, help='Path to output BAM file')
    args = parser.parse_args()

    fix_primary_alignments(args.input, args.output)


if __name__ == '__main__':
    main()
