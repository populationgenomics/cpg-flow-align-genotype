import numpy as np
import pandas as pd
from cpg_utils import to_path

NEGATIVES = ['Not flagged', 'Not applicable', 'None', None, pd.NA, np.nan]


def scan_vntyper_html(html_file: str) -> dict[str, bool]:
    """
    Parse a VNtyper HTML file to detect positive/negative calls.
    Positive status is implied by the relevant table containing a non-negative row.
    """

    with to_path(html_file).open() as f:
        df_list = pd.read_html(f)

    kestrel_positive = bool(np.any(~(df_list[1]['Variant'].isin(NEGATIVES))))
    advntr_positive = bool(np.any(~(df_list[2]['Variant'].isin(NEGATIVES))))

    return {'advntr': advntr_positive, 'kestrel': kestrel_positive}
