"""
This script is set to run at the end of the align/genotype workflow, and only runs for projects which are configured
to run VNtyper

It queries Metamist for all analysis entries representing VNtyper results, reduces those down to 1/sequencing group

This will run scoped by dataset and sequencing_type.
"""

from argparse import ArgumentParser
from pathlib import Path

import jinja2

from cpg_utils import config, to_path
from metamist import graphql


JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'

REPORT_QUERY = graphql.gql(
    """
    query ReportQuery($project: String!, $metaFilter: JSON) {
      project(name: $project) {
        analyses(type: {eq: "web"}, meta: $metaFilter) {
          id
          meta
          output
          sequencingGroups {
            id
          }
          timestampCompleted
        }
      }
    }"""
)
REPORT_EMOJI: dict[bool, str] = {True: '🟢', False: '🔴'}
WEB_BASE: str = 'gs://cpg-{}-main-web'
WEB_URL_BASE: str = 'https://main-web.populationgenomics.org.au/{}'


def translate_url(og_url: str) -> str:
    """Translate the absolute location of the HTML document to a proxy-rendered link."""


def query_for_reports(dataset: str, sequencing_type: str) -> dict[str, dict[str, str]]:
    """
    Execute a GQL query for all relevant VNtyper reports. Minimise to one-per-SGID.

    Args:
        dataset: str, the name of the dataset/project
        sequencing_type: str, the type of assay (typicall exome/genome)

    Returns:
        dict of SGID: Report
    """

    report_lookup: dict[str, dict[str, bool | str]] = {}

    # build a meta query object
    meta_param = {
        'sequencing_type': sequencing_type,
        'stage': 'RunVntyper',
    }

    results = graphql.query(REPORT_QUERY, variables={'project': dataset, 'meta': meta_param})

    # trim `-test` off the dataset string, to present `-test-test` in the URL
    dataset_string = dataset.replace('-test', '')

    access_level = config.config_retrieve(['workflow', 'access_level'])

    # populate the expected URL portion, and adjust if test
    web_base = WEB_BASE.format(dataset_string)
    if access_level == 'test':
        web_base = web_base.replace('main', 'test')

    # populate the proxy-enabled URL portion, and adjust if test
    proxy_base = WEB_URL_BASE.format(dataset_string)
    if access_level == 'test':
        proxy_base = proxy_base.replace('main', 'test')

    # iterate over the reports and grab them all
    for report in results['project']['analyses']:
        meta_dict = report['meta']

        # pull out the URL from the analysis entry, swap the real URL for a proxy URL
        url = report['output'].replace(web_base.format(dataset), proxy_base.format(dataset))

        # capture a minrep of all the fields we want to present in the report
        report_lookup[report['sequencingGroups'][0]['id']] = {
            'advntr': REPORT_EMOJI[meta_dict['advntr']],
            'kestrel': REPORT_EMOJI[meta_dict['kestrel']],
            'participant': meta_dict['participant_id'],
            'timestamp': report['timestampCompleted'].split('T')[0],
            'url': url,
        }

    return report_lookup


def main(dataset: str, output: str) -> None:
    sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
    template_context = {
        'title': f'VNtyper index for {dataset}, {sequencing_type}',
        'reports': query_for_reports(dataset, sequencing_type),
    }

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('vntyper_index.html.jinja')
    content = template.render(**template_context)

    # write to common web bucket - either attached to a single dataset, or communal
    print(f'Writing {template_context["title"]} to {output}')
    to_path(output).write_text('\n'.join(line for line in content.split('\n') if line.strip()))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True, help='Dataset to generate index page for.')
    parser.add_argument('--output', required=True, help='Path to write new HTML file to.')
    args = parser.parse_args()
    main(args.dataset, args.output)
