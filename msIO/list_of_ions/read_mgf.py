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

import logging

logger = logging.getLogger(__name__)


def ignore_line(line: str) -> bool:
    # skip comments and blank lines
    return (len(line) == 0) or any([line.startswith(x) for x in '#;!/'])


class MgfImportManager(BaseLib, FeatureManager):
    """create dict of features to spectra by parsing an mgf file"""

    @property
    def feature_ids(self) -> np.ndarray:
        return self._feature_ids

    def __init__(self, path_mgf: str, accumulate_spectra: bool = True, mz_resolution: int = None):
        """
        accumulate_spectra: if this is set to True, will merge spectra to the provided resolution
        mz_resolution: mz/dmz up to which peaks are distinguishable (used for accumulating spectra)
        filter_ms1_spectra: if this is set to True,
        """

        self.peak_list: list[PeakList] = []
        _feature_ids: list[int] = []
        _ms_level: list[int] = []
        _ions: list[str] = []

        header = []
        entries: list[dict] = []
        is_ion: bool = False  # could contain header

        # get line count
        with open(path_mgf, 'rb') as f:
            n_lines = sum(1 for _ in f)

        with open(path_mgf, 'r') as f:
            lines_ion = []
            for line in tqdm(f, desc='reading mgf file', total=n_lines):
                if line.startswith('BEGIN IONS'):
                    is_ion = True
                elif line.startswith('END IONS'):
                    feature_props: dict = parse_ion_props(lines_ion)
                    peak_list = PeakList.from_lines(lines_ion)
                    self.peak_list.append(peak_list)
                    entries.append(feature_props)
                    _feature_ids.append(int(feature_props['feature_id']))
                    _ms_level.append(int(feature_props['ms_level']))
                    _ions.append(feature_props.get('ion', '[M+?]'))
                    lines_ion = []
                    is_ion = False
                elif is_ion:
                    lines_ion.append(line)
                elif not ignore_line(line):
                    header.append(line)

        self.df_features: pd.DataFrame = pd.DataFrame(entries)
        # self.df_features.set_index('feature_id')
        if (('rt_minutes' not in self.df_features.columns)
                and ('rt_seconds' in self.df_features.columns)):
            self.df_features.loc[:, 'rt_minutes'] = self.df_features.rt_seconds / 60
        elif (('rt_seconds' not in self.df_features.columns)
              and ('rt_minutes' in self.df_features.columns)):
            self.df_features.loc[:, 'rt_seconds'] = self.df_features.rt_minutes * 60
        if 'polarity' not in self.df_features.columns:
            assert 'charge' in self.df_features.columns, 'charge or polarity attribute is required'
            self.df_features.loc[:, 'polarity'] = self.df_features.charge.apply(lambda x: "POSITIVE" if x > 0 else "NEGATIVE")

        if accumulate_spectra:
            assert mz_resolution is not None, 'mz_resolution must be provided if accumulate_spectra is True'
            # merge spectra with the same feature id, ms_level and ion
            ms_spec_keys = [tuple(row) for _, row in self.df_features.loc[:, ['feature_id', 'ms_level', 'ion']].iterrows()]
            unique_ms_spec_keys = set(ms_spec_keys)

            peak_lists_merged: list[PeakList] = []
            for ms_spec_key in tqdm(unique_ms_spec_keys, desc='merging spectra', total=len(unique_ms_spec_keys)):
                peak_lists = [pl for k, pl in zip(ms_spec_keys, self.peak_list) if k == ms_spec_key]
                # merge
                pl_merged: PeakList = peak_lists[0]
                for pl in peak_lists[1:]:
                    pl_merged = pl_merged.merge_with(pl, mz_resolution=mz_resolution)
                peak_lists_merged.append(pl_merged)
            self.peak_list = peak_lists_merged
            # some columns will no longer have the correct values but they are not relevant for the MgfFeatures anyway
            self.df_features.drop_duplicates(subset=['feature_id', 'ms_level', 'ion'], inplace=True)

            self._peak_dict: dict[tuple[int, int, str], PeakList] = dict(zip(
                unique_ms_spec_keys,
                self.peak_list)
            )
        else:
            if self.df_features.loc[:, ['feature_id', 'ms_level', 'ion']].duplicated().any():
                logger.warning(
                    'Found more than one MS spectrum per level for some features, consider setting '
                    'accumulate_spectra=True'
                )

            self._peak_dict: dict[tuple[int, int, str], PeakList] = dict(zip(
                zip(_feature_ids, _ms_level, _ions),
                self.peak_list)
            )
        self._feature_ids: np.ndarray[int] = np.unique(_feature_ids)

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
            props.pop('polarity', None)
            props.pop('feature_id')
            ms_specs.append(MsSpec(peaks=peaks, **{k: v for k, v in props.items() if k in MsSpec.__annotations__}))

        f = FeatureMgf(
            feature_id=f_id,
            polarity=df_sub.polarity.iat[0],
            ms_specs=ms_specs
        )
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
    # path_mgf_sirius = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.sirius.mgf"
    path_mgf = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\mzmine\guaymas_new_sirius.mgf"

    unmerged_mgf = MgfImportManager(path_mgf, accumulate_spectra=False)
    mgf = MgfImportManager(path_mgf, accumulate_spectra=True, mz_resolution=40_000)

    # peak_lists = mgf.get_ms2(mz=636.53379)
    #
    # peak_lists[1][1].plot()
    # plt.show()

    fig, axs = plt.subplots(ncols=2, sharex=True, sharey=True)
    mgf._peak_dict[(35, 2, np.nan)].plot(ax=axs[0])
    for i, idx in enumerate([5, 6, 7, 8, 9, 10]):
        unmerged_mgf.peak_list[idx].plot(ax=axs[1], linefmt=f'C{i}')

    f = mgf.get_feature(1)
