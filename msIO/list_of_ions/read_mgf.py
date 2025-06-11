"""This module allows access to exported MS2 spectra from mgf files"""
from dataclasses import dataclass
from typing import Literal, Self

import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

from .base import BaseLib, PeakList


def _parse_ion_props(inpt: list[str]) -> dict:
    kwargs: dict[str, int | float | str] = {}
    for line in inpt:
        if line.startswith('FEATURE_ID'):
            kwargs['feature_id'] = int(line.split('=')[1])
        elif line.startswith('PEPMASS'):
            kwargs['mz'] = float(line.split('=')[1])
        elif line.startswith('MSLEVEL'):
            kwargs['ms_level'] = int(line.split('=')[1])
        elif line.startswith('CHARGE'):
            # sign is trailing
            charge = line.split('=')[1].strip('\n')
            sign = charge[-1]
            mag = charge[:-1]
            kwargs['charge'] = int(sign + mag)
        elif line.startswith('POLARITY'):
            kwargs['polarity'] = line.split('=')[1].strip('\n').lower()
        elif line.startswith('ION'):
            kwargs['ion'] = line.split('=')[1].strip('\n')
        elif line.startswith('RTINMINUTES'):
            kwargs['rt_minutes'] = float(line.split('=')[1])
        elif line.startswith('RTINSECONDS'):
            kwargs['rt_seconds'] = float(line.split('=')[1])
        elif line.startswith('TITLE'):
            kwargs['title'] = line.split('=',1)[1]
        else:
            continue
    return kwargs


@dataclass
class Entry:
    feature_id: int
    mz: float
    ms_level: int
    charge: int
    polarity: Literal['pos', 'neg']
    ion: str
    rt_seconds: float = None
    rt_minutes: float = None

    def __post_init__(self):
        if not self.rt_seconds and self.rt_minutes:
            self.rt_seconds = self.rt_minutes * 60

        if not self.rt_minutes and self.rt_seconds:
            self.rt_minutes = self.rt_seconds / 60

    @classmethod
    def from_lines(cls, inpt: list[str]) -> Self:
        return cls(**_parse_ion_props(inpt))


class MGFReader(BaseLib):
    """create dict of features to spectra by parsing an mgf file"""
    def __init__(self, path_mgf: str):
        self.peak_list: list[PeakList] = []
        entries: list[dict] = []
        with open(path_mgf, 'r') as f:
            lines_ion = []
            for line in tqdm(f, desc='reading mgf file'):
                lines_ion.append(line)
                if line.startswith('END IONS'):
                    feature_props: dict = _parse_ion_props(lines_ion)
                    peaks = PeakList.from_lines(lines_ion)
                    self.peak_list.append(peaks)
                    entries.append(feature_props)
                    lines_ion = []

        self.features: pd.DataFrame = pd.DataFrame(entries)
        # self.features.set_index('feature_id')
        if ('rt_minutes' not in self.features.columns) and ('rt_seconds' in self.features.columns):
            self.features.loc[:, 'rt_minutes'] = self.features.rt_seconds / 60
        elif ('rt_seconds' not in self.features.columns) and ('rt_minutes' in self.features.columns):
            self.features.loc[:, 'rt_seconds'] = self.features.rt_minutes * 60

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

    mgf = MGFReader(path_mgf_sirius)

    peak_lists = mgf.get_ms2(mz=636.53379)

    peak_lists[1][1].plot()
    plt.show()
