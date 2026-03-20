import numpy as np
import pandas as pd

from msIO.feature_managers.base import FeatureManager
from msIO.features.metaboscape import FeatureMetaboScape, METABOSCAPE_CSV_RENAME_COLUMNS
from msIO.environmental.sample import Sample


class MetaboscapeImportManager(FeatureManager):
    def __init__(self, path_metaboscape_export_file: str, path_metaboscape_clipboard_file: str = None):
        self.path_metaboscape_export_file = path_metaboscape_export_file

        _df = pd.read_csv(path_metaboscape_export_file).rename(columns=METABOSCAPE_CSV_RENAME_COLUMNS)
        # discard mean, max intensity columns
        drop_cols = []
        for col in _df.columns:
            if col.endswith('Intensity'):
                drop_cols.append(col)
        self._df: pd.DataFrame = _df.drop(columns=drop_cols)

        if path_metaboscape_clipboard_file is not None:
            self._add_from_clipboard(path_metaboscape_clipboard_file)

        self._sample_name_to_sample: dict[str, Sample] = {}

    def _add_from_clipboard(self, path_file: str):
        # rows should match 1 to 1
        df_clip = pd.read_csv(path_file, sep='\t')

        err_msg = 'data in rows of the export and clipboard files should be in the same order'
        assert df_clip.shape[0] == self._df.shape[0], err_msg
        # check the order
        assert np.all(np.abs(self._df.loc[:, 'rt_seconds'] / 60 - df_clip.loc[:, 'RT [min]']) < .02), err_msg

        columns_to_transfer = ['m/z meas.', 'Flags', 'AQ', 'Annotation Source', 'Include']

        for c in columns_to_transfer:
            self._df.loc[:, c] = df_clip.loc[:, c]
        self._df.rename(columns=METABOSCAPE_CSV_RENAME_COLUMNS, inplace=True)

    @property
    def feature_ids(self):
        return self._df.feature_id.unique()

    def _inner_missing_feature(self, f_id) -> None:
        idx = np.argwhere(self._df.feature_id == f_id)[0][0]
        f = FeatureMetaboScape.from_dataframe_row(self._df.iloc[idx, :], self._sample_name_to_sample)
        self._features[f_id] = f


if __name__ == '__main__':
    # path_file_csv = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.csv"
    path_file_csv = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas really new method\MetaboScape\and a backup.csv"
    path_file_clipboard = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas really new method\MetaboScape\from_clipboard.csv"

    mm = MetaboscapeImportManager(path_file_csv, path_file_clipboard)

    f = mm.get_feature(397)
