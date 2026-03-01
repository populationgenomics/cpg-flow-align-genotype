"""
Script which reads a .interval_list file and a fasta.dict
and rewrites the .interval_list file with the hashes from the fasta.dict file.
"""

import argparse
import re

from cpg_utils import to_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Rewrite the sequence dictionary in an interval list file using the hashes from a fasta.dict file."""
    )
    parser.add_argument('interval_list', help='Path to the input .interval_list file.')
    parser.add_argument('fasta_dict', help='Path to the input fasta.dict file.')
    parser.add_argument('output_interval_list', help='Path to the output .interval_list file.')
    return parser.parse_args()


def read_fasta_dict(fasta_dict_path):
    """Read the fasta.dict file and return a dictionary mapping sequence names to their hashes."""
    fasta_dict = {}
    with open(to_path(fasta_dict_path)) as f:
        for line in f:
            if not line.startswith('@SQ'):
                continue
            # Lines look like: "@SQ	SN:chr1	LN:XXXX	M5:<hash>	UR:file:/path/to/ref.fasta"
            seq_name, seq_hash = re.search(r'SN:(\S+)', line).group(1), re.search(r'M5:(\S+)', line).group(1)
            fasta_dict[seq_name] = seq_hash
    return fasta_dict


def rewrite_interval_list(interval_list_path, fasta_dict, output_path):
    with open(to_path(interval_list_path)) as f:
        lines = f.readlines()
    with to_path(output_path).open('w') as f:
        for line in lines[:]:
            if line.startswith('@SQ'):
                seq_name = re.search(r'SN:(\S+)', line).group(1)
                current_hash = re.search(r'M5:(\S+)', line).group(1)
                if seq_name in fasta_dict:
                    new_hash = fasta_dict[seq_name]
                    if new_hash != current_hash:
                        new_line = line.replace(f'M5:{current_hash}', f'M5:{new_hash}')
                        print(f'Updated hash for {seq_name}: {current_hash} -> {new_hash}')
                        f.write(new_line)
                        continue

            f.write(line)


if __name__ == '__main__':
    args = parse_args()
    fasta_dict = read_fasta_dict(args.fasta_dict)
    rewrite_interval_list(args.interval_list, fasta_dict, args.output_interval_list)
    print('Done.')
