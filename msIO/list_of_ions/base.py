import warnings
from typing import Self, Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class PeakList:
    """parse list of ions to create spectrum"""

    def __init__(
            self,
            mzs: Iterable[float] = None,
            intensities: Iterable[float] = None,
            peaks: dict[float, float] = None
    ) -> None:
        if (mzs is None) and (intensities is None) and (peaks is None):
            self.peaks = {}
        elif peaks is None:
            assert (mzs is not None) and (intensities is not None), 'provide either peaks OR (mzs AND intensities)'
            assert len(mzs) == len(intensities)
            self.peaks = {}
            for mz, i in zip(mzs, intensities):
                self.peaks[mz] = i
        else:
            assert (mzs is None) and (intensities is None), 'provide either peaks OR (mzs AND intensities)'
            self.peaks = peaks

    def __add__(self, other: Self) -> Self:
        new_peaks = self.peaks.copy()
        for mz, i in other.peaks.items():
            if mz not in new_peaks:
                new_peaks[mz] = 0
            new_peaks[mz] += i
        return self.__class__(peaks=new_peaks)

    @classmethod
    def from_lines(cls, inpt: list[str], splitter=' '):
        lines_peaks = [l for l in inpt if l[0].isnumeric()]

        mzs = []
        ints = []
        for line in lines_peaks:
            mz, i = line.split(splitter)
            mzs.append(float(mz))
            ints.append(float(i))

        return cls(mzs, ints)

    def plot(self, ax: plt.Axes = None) -> plt.Axes:
        if ax is None:
            _, ax = plt.subplots()

        ax.stem(self.peaks.keys(), self.peaks.values(), markerfmt='')
        ax.set_xlabel('m/z in Da')
        ax.set_ylabel('Intensity')

        return ax


class BaseLib:
    df_features: pd.DataFrame = None
    peak_list: list[PeakList] = None

    def _get_ms2(
            self,
            mz: float = None,
            rt_minutes: float = None,
            rt_seconds: float = None,
            mass_tolerance: float = None,
            rt_minutes_tolerance: float = None,
            rt_seconds_tolerance: float = None,
            ion: str = None,
            feature_index: int = None
    ) -> tuple[pd.DataFrame, list[PeakList]]:
        if self.peak_list is None:
            raise AttributeError('peak_list must be initialized first')

        """Fetch the MS2 spectrum for a specific mz and RT"""
        assert (rt_minutes is None) or (rt_seconds is None), \
            'give RT either in seconds or minutes, but not both'

        if mz is not None:
            assert mass_tolerance is not None, \
                'mass tolerance is required if pepmass is provided'
        if rt_minutes is not None:
            assert rt_minutes_tolerance is not None, \
                'rt_minutes_tolerance is required if rt_minutes is provided'
        if rt_seconds is not None:
            assert rt_seconds_tolerance is not None, \
                'rt_seconds_tolerance is required if rt_seconds is provided'

        mask = self.df_features.ms_level == 2
        if mz is not None:
            mask &= np.abs(self.df_features.mz - mz) < mass_tolerance
        if rt_minutes is not None:
            mask &= np.abs(self.df_features.rt_minutes - rt_minutes) < .01
        if rt_seconds is not None:
            mask &= np.abs(self.df_features.rt_seconds - rt_seconds) < .002

        # pick right adduct
        if ion is not None:
            mask = (self.df_features.ion == ion) & mask
        else:  # prefer H, then NH4, then Na adduct
            for ion in ('[M+H]+', '[M+NH4]+', '[M+Na]+'):
                if ion in self.df_features.ion[mask]:
                    mask = (self.df_features.ion == ion) & mask
                    break
            else:
                warnings.warn('unable to find any of the standard adducts, '
                              'using any adduct that is available')

        n_matches = mask.sum()
        if n_matches != 1:
            if rt_minutes is not None:
                rt = f'{rt_minutes:.2f} min'
            elif rt_seconds is not None:
                rt = f'{rt_seconds:.1f} sec'
            else:
                rt = ''
            if rt != '':
                rt = 'and ' + rt
            warnings.warn(f'found {n_matches} matches for {mz} Da {rt}')

        ids = self.df_features.index[mask].to_list()
        peak_lists: list[PeakList] = []
        for _id in ids:
            peak_lists.append(self.peak_list[_id])
        return self.df_features.loc[mask, :], peak_lists
