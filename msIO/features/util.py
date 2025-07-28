from typing import Callable

import pandas as pd

from msIO.features.base import FeatureBaseClass


def read_csv(file_csv: str, column_renamer: dict[str, str], hook: Callable = None, **kwargs):
    """Wrapper for pandas read_csv"""
    if hook is None:
        hook = lambda x: x

    df = pd.read_csv(file_csv, **kwargs)
    df.rename(columns=column_renamer, inplace=True)
    return hook(df)


def list_of_features_to_dataframe(features: list[FeatureBaseClass]) -> pd.DataFrame:
    columns = set()
    for attrs in features.__dict__.keys():
        columns |= set(attrs)
    ...
