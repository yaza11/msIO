import warnings
from typing import Self, Iterable, Optional, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sqlalchemy import ForeignKey, String, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship
from enum import Enum as PyEnum
from sqlalchemy import Enum

from msIO.features.base import SqlBaseClass, FeatureBaseClass


class SpectrumType(PyEnum):
    retention_time = "retention_time"
    mobility = "mobility"
    mass_over_charge = "mass_over_charge"
    mass_over_charge_isotopes = "mass_over_charge_isotopes"
    mass_over_charge_fragments = "mass_over_charge_fragments"


class PeakFeature(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "peak"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)

    peak_list_id: Mapped[int] = mapped_column(ForeignKey("peak_list.id"))
    peak_list: Mapped["PeakList"] = relationship(back_populates="peaks")  # every peak is part of a peak list


class PeakList(SqlBaseClass, FeatureBaseClass):
    """parse list of ions to create spectrum"""
    __tablename__ = "peak_list"

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    spectrum_type: Mapped[SpectrumType | None] = mapped_column(
        Enum(SpectrumType, name='spectrum_type_enum'),
        nullable=True
    )

    peaks: Mapped[list["PeakFeature"]] = relationship(
        back_populates="peak_list",
        cascade="all, delete-orphan"
    )

    def __init__(
            self,
            mzs: Optional[Iterable[float]] = None,
            intensities: Optional[Iterable[float]] = None,
            peaks: Optional[Union[dict[float, float], Iterable[PeakFeature]]] = None,
            name: Optional[str] = None
    ) -> None:
        self.name = name
        self.peaks = self._build_peaks(mzs, intensities, peaks)

    def _build_peaks(
            self,
            mzs: Optional[Iterable[float]] = None,
            intensities: Optional[Iterable[float]] = None,
            peaks: Optional[Union[dict[float, float], Iterable[PeakFeature]]] = None
    ) -> list[PeakFeature]:
        if (mzs is None) and (intensities is None) and (peaks is None):
            return []

        if peaks is None:
            assert (mzs is not None) and (intensities is not None), 'provide either peaks OR (mzs AND intensities)'
            assert len(mzs) == len(intensities)
            peaks = []
            for mz, i in zip(mzs, intensities):
                peaks.append(PeakFeature(mz=mz, intensity=i))
            return [PeakFeature(mz=p.mz, intensity=p.intensity, peak_list=self) for p in peaks]

        assert (mzs is None) and (intensities is None), 'provide either peaks OR (mzs AND intensities)'
        if isinstance(peaks, dict):
            return [PeakFeature(mz=k, intensity=v, peak_list=self) for k, v in peaks.items()]
        return [p
                if p.peak_list is self
                else PeakFeature(mz=p.mz, intensity=p.intensity, peak_list=self)
                for p in peaks]

    @property
    def mzs(self) -> list[float]:
        return [p.mz for p in self.peaks]

    @property
    def intensities(self) -> list[float]:
        return [p.intensity for p in self.peaks]

    def __add__(self, other: Self) -> Self:
        new_peaks = self.peaks.copy()
        mzs = self.mzs
        for p in other.peaks:
            mz = p.mz
            i = p.intensity
            if mz not in mzs:
                new_peaks.append(PeakFeature(mz=mz, intensity=0.))
            idx = mzs.index(mz)
            new_peaks[idx].intensity += i
        return self.__class__(peaks=new_peaks)

    @classmethod
    def from_lines(cls, inpt: list[str], splitter=' ', name: Optional[str] = None) -> Self:
        lines_peaks = [l for l in inpt if l[0].isnumeric()]

        mzs = []
        ints = []
        for line in lines_peaks:
            mz, i = line.split(splitter)[:2]
            mzs.append(float(mz))
            ints.append(float(i))

        return cls(mzs, ints, name=name)

    def plot(self, ax: plt.Axes = None) -> plt.Axes:
        if ax is None:
            _, ax = plt.subplots()

        ax.stem(self.mzs, self.intensities, markerfmt='')
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


if __name__ == '__main__':
    peak = PeakFeature(mz=100, intensity=10)

    pl = PeakList(mzs=[1, 2, 3], intensities=[3, 4, 5])
