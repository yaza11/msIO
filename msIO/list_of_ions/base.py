import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



class PeakList:
    """parse list of ions to create spectrum"""
    def __init__(self, inpt: list[str], splitter=' '):
        self.peaks: dict[float, float] = {}
        for line in inpt:
            mz, ints = line.split(splitter)
            self.peaks[float(mz)] = float(ints)

    @classmethod
    def from_lines(cls, inpt: list[str], **kwargs):
        lines_peaks = [l for l in inpt if l[0].isnumeric()]
        return cls(lines_peaks, **kwargs)

    def plot(self, ax: plt.Axes = None) -> plt.Axes:
        if ax is None:
            _, ax = plt.subplots()

        ax.stem(self.peaks.keys(), self.peaks.values(), markerfmt='')
        ax.set_xlabel('m/z in Da')
        ax.set_ylabel('Intensity')

        return ax


class BaseLib:
    features: pd.DataFrame = None
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

        mask = self.features.ms_level == 2
        if mz is not None:
            mask &= np.abs(self.features.mz - mz) < mass_tolerance
        if rt_minutes is not None:
            mask &= np.abs(self.features.rt_minutes - rt_minutes) < .01
        if rt_seconds is not None:
            mask &= np.abs(self.features.rt_seconds - rt_seconds) < .002

        # pick right adduct
        if ion is not None:
            mask = (self.features.ion == ion) & mask
        else:  # prefer H, then NH4, then Na adduct
            for ion in ('[M+H]+', '[M+NH4]+', '[M+Na]+'):
                if ion in self.features.ion[mask]:
                    mask = (self.features.ion == ion) & mask
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

        ids = self.features.index[mask].to_list()
        peak_lists: list[PeakList] = []
        for _id in ids:
            peak_lists.append(self.peak_list[_id])
        return self.features.loc[mask, :], peak_lists
