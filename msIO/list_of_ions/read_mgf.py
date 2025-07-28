"""
Read mgf (Mascot Generic Format) files

https://www.matrixscience.com/help/data_file_help.html
"""
from dataclasses import dataclass
from typing import Literal, Self, Iterable

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

from msIO.feature_managers.base import FeatureManager
from msIO.features.mgf import parse_ion_props, FeatureMgf, MsSpec
from msIO.list_of_ions.base import BaseLib, PeakList


def ignore_line(line: str) -> bool:
    # skip comments and blank lines
    return (len(line) == 0) or any([line.startswith(x) for x in '#;!/'])


class MgfImportManager(BaseLib, FeatureManager):
    """create dict of features to spectra by parsing an mgf file"""

    @property
    def feature_ids(self) -> np.ndarray:
        return self._feature_ids

    def __init__(self, path_mgf: str):
        self.peak_list: list[PeakList] = []
        _feature_ids: list[int] = []
        _ms_level: list[int] = []
        _ions: list[str] = []

        header = []
        entries: list[dict] = []
        is_ion: bool = False  # could contain header
        with open(path_mgf, 'r') as f:
            lines_ion = []
            for line in tqdm(f, desc='reading mgf file'):
                if line.startswith('BEGIN IONS'):
                    is_ion = True
                elif line.startswith('END IONS'):
                    feature_props: dict = parse_ion_props(lines_ion)
                    peaks = PeakList.from_lines(lines_ion)
                    self.peak_list.append(peaks)
                    entries.append(feature_props)
                    _feature_ids.append(int(feature_props['feature_id']))
                    _ms_level.append(int(feature_props['ms_level']))
                    _ions.append(feature_props['ion'])
                    lines_ion = []
                    is_ion = False
                elif is_ion:
                    lines_ion.append(line)
                elif not ignore_line(line):
                    header.append(line)

        self.df_features: pd.DataFrame = pd.DataFrame(entries)
        # self.df_features.set_index('feature_id')
        if ('rt_minutes' not in self.df_features.columns) and ('rt_seconds' in self.df_features.columns):
            self.df_features.loc[:, 'rt_minutes'] = self.df_features.rt_seconds / 60
        elif ('rt_seconds' not in self.df_features.columns) and ('rt_minutes' in self.df_features.columns):
            self.df_features.loc[:, 'rt_seconds'] = self.df_features.rt_minutes * 60

        self._feature_ids: np.ndarray[int] = np.unique(_feature_ids)
        self._peak_dict: dict[tuple[int, int, str], PeakList] = dict(zip(
            zip(_feature_ids, _ms_level, _ions),
            self.peak_list)
        )

    def _inner_missing_feature(self, f_id) -> None:
        # get properties from dataframe
        mask_id = self.df_features.feature_id == f_id
        if mask_id.sum() == 0:
            return
        df_sub = self.df_features.loc[mask_id, :]

        ms_specs: list[MsSpec] = []
        for idx, row in df_sub.iterrows():
            # create keys to check for which ones we have MS spectra
            key = row.feature_id, row.ms_level, row.ion
            if key not in self._peak_dict:
                continue
            peaks = self._peak_dict[key]
            props = row.to_dict()
            props.pop('polarity')
            props.pop('feature_id')
            ms_specs.append(MsSpec(peaks=peaks, **props))

        f = FeatureMgf(feature_id=f_id, polarity=df_sub.polarity.iat[0], ms_specs=ms_specs)
        self._features[f_id] = f

    def get_ms2(
            self,
            mz: float = None,
            rt_minutes: float = None,
            rt_seconds: float = None,
            mass_tolerance: float = 1e-3,
            rt_minutes_tolerance: float = .01,
            rt_seconds_tolerance: float = .002
    ) -> tuple[pd.DataFrame, list[PeakList]]:
        return self._get_ms2(
            mz, rt_minutes, rt_seconds, mass_tolerance,
            rt_minutes_tolerance, rt_seconds_tolerance
        )


if __name__ == '__main__':
    path_mgf_sirius = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.sirius.mgf"

    mgf = MgfImportManager(path_mgf_sirius)

    # peak_lists = mgf.get_ms2(mz=636.53379)
    #
    # peak_lists[1][1].plot()
    # plt.show()

    f = mgf.get_feature(2)
