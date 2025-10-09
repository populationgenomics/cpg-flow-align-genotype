#!/usr/bin/env python3

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow

from align_genotype.dataset_stages import CramMultiQC, GvcfMultiQC, SomalierPedigree


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    stages = [SomalierPedigree, CramMultiQC, GvcfMultiQC]

    run_workflow(name='align_genotype', stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
