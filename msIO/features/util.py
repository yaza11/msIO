from typing import Callable

import pandas as pd


def read_csv(file_csv: str, column_renamer: dict[str, str], hook: Callable = None, **kwargs):
    """Wrapper for pandas read_csv"""
    if hook is None:
        hook = lambda x: x

    df = pd.read_csv(file_csv, **kwargs)
    df.rename(columns=column_renamer, inplace=True)
    return hook(df)
