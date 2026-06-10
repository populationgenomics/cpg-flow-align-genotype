import pandas as pd

from cpg_utils import to_path

NEGATIVES = ['Not flagged', 'Not applicable', 'None']


def scan_vntyper_html(html_file: str) -> dict[str, bool]:
    """
    Parse a VNtyper HTML file to detect positive/negative calls.
    Positive status is implied by the relevant table containing a non-negative row.
    """

    with to_path(html_file).open() as f:
        df_list = pd.read_html(f)

    advntr_positive = False
    kestrel_positive = False

    # kestrel
    for _, row in df_list[1].iterrows():
        if row['Variant'] not in NEGATIVES:
            kestrel_positive = True
            break

    # advntr
    for _, row in df_list[2].iterrows():
        if row['Variant'] not in NEGATIVES:
            advntr_positive = True
            break

    return {'advntr': advntr_positive, 'kestrel': kestrel_positive}
