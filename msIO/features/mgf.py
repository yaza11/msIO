from dataclasses import dataclass
from typing import Literal, Self

import numpy as np

from msIO import PeakList
from msIO.features.base import FeatureBaseClass


# select which ion is most likely (first in list)
ION_PREFERENCES = ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+', '[M]+', '[M+H+H]2+', '[M+H+H2]3+']


def parse_ion_props(inpt: list[str]) -> dict:
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
            kwargs['title'] = line.split('=', 1)[1]
        else:
            continue
    return kwargs


@dataclass
class MsSpec:
    mz: float = None
    ms_level: int = None
    charge: int = None

    rt_seconds: float = None
    ion: str = None
    rt_minutes: float = None

    peaks: PeakList = None


@dataclass
class FeatureMgf(FeatureBaseClass):
    """Multiple entries can belong to a single feature id
    (when deconvolution fails or when multiple adducts are found).
    This object sets ambiguous parameters by choosing from the preferred ion
    order"""
    feature_id: int = None
    polarity: Literal['pos', 'neg'] = None  # should always be the same for all adducts
    ms_specs: list[MsSpec] = None
    has_multiple_adducts: bool = None

    # set from preferred ion
    mz: float = None
    charge: int = None
    rt_seconds: float = None
    ion = None
    rt_minutes: float = None

    ms1: PeakList = None
    ms2: PeakList = None

    @classmethod
    def from_lines(cls, inpt: list[str]) -> Self:
        return cls(**parse_ion_props(inpt))

    def __post_init__(self) -> None:
        """Set properties from preferred ion"""
        if (self.ms_specs is None) or len(self.ms_specs) == 0:
            self.has_multiple_adducts = False
            return

        adducts: dict[tuple[str, int], MsSpec] = {(msspec.ion, msspec.ms_level): msspec for msspec in self.ms_specs}
        self.has_multiple_adducts = len(set([k[0] for k in adducts.keys()])) > 1
        # on first pass, try to find adduct for which both levels exist
        for add_pref in ION_PREFERENCES:
            if ((add_pref, 2) in adducts) and ((add_pref, 1) in adducts):
                key = (add_pref, 2)  # prefer MS2
                ms_pref = adducts[key]
                break
        else:  # try to find just 1 level
            for add_pref in ION_PREFERENCES:
                # prefer MS2
                if (key := (add_pref, 2)) in adducts:
                    ms_pref = adducts[key]
                    break
                elif (key := (add_pref, 1)) in adducts:
                    ms_pref = adducts[key]
                    break
            else:
                # pick any of those that are there
                key = list(adducts.keys())[0]
                add_pref = key[0]
                ms_pref = adducts[key]

        update_props = ['ion', 'rt_seconds', 'rt_minutes', 'charge', 'mz']
        self.__dict__ |= {k: v
                          for k, v in ms_pref.__dict__.items()
                          if k in update_props}

        for level in [1, 2]:
            if (k := (add_pref, level)) in adducts:
                self.__setattr__(f'ms{level}', adducts[k].peaks)


def test():
    path_mgf_sirius = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.sirius.mgf"

    lines_ion = []
    is_ion = False
    with open(path_mgf_sirius, 'r') as f:
        for line in f:
            if line.startswith('BEGIN IONS'):
                is_ion = True
            elif line.startswith('END IONS'):
                break
            elif is_ion:
                lines_ion.append(line)

    f = FeatureMgf.from_lines(lines_ion)
    return f


if __name__ == '__main__':
    f = test()
