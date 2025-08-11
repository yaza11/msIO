from typing import Iterable

import numpy as np
import pandas as pd
from tqdm import tqdm

from msIO.feature_managers.base import FeatureManager
from msIO.features.metaboscape import FeatureMetaboScape, METABOSCAPE_CSV_RENAME_COLUMNS


class MetaboscapeImportManager(FeatureManager):
    def __init__(self, path_metaboscape_export_file: str):
        self.path_metaboscape_export_file = path_metaboscape_export_file

        _df = pd.read_csv(path_metaboscape_export_file).rename(columns=METABOSCAPE_CSV_RENAME_COLUMNS)
        # discard mean, max intensity columns
        drop_cols = []
        for col in _df.columns:
            if col.endswith('Intensity'):
                drop_cols.append(col)
        self._df: pd.DataFrame = _df.drop(columns=drop_cols)

    @property
    def feature_ids(self):
        return self._df.feature_id.unique()

    def _inner_missing_feature(self, f_id) -> None:
        idx = np.argwhere(self._df.feature_id == f_id)[0][0]
        f = FeatureMetaboScape.from_dataframe_row(self._df.iloc[idx, :])
        self._features[f_id] = f


if __name__ == '__main__':
    path_file_csv = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.csv"

    mm = MetaboscapeImportManager(path_file_csv)

    f = mm.get_feature(1)
