from dataclasses import dataclass
from typing import Self, Optional

from enum import Enum as PyEnum
from sqlalchemy.orm import Mapped, mapped_column, relationship
from sqlalchemy import ForeignKey, Boolean, String, Integer, Float
from sqlalchemy import Enum

from msIO.features.base import FeatureBaseClass, SqlBaseClass
from msIO.list_of_ions.base import PeakList

# select which ion is most likely (first in list)
ION_PREFERENCES = ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+', '[M]+', '[M+H+H]2+', '[M+H+H2]3+']


def parse_ion_props(inpt: list[str]) -> dict:
    """Parse mgf entry"""
    do_nothing = lambda x: x
    keyword_to_label_and_conversion = {
        'feature_id': ('feature_id', int),
        'pepmass': ('mz', float),
        'mslevel': ('ms_level', int),
        'charge': ('charge', lambda c: int(c[-1] + c[:-1])),  # e.g. 2- --> -2
        'rtinminutes': ('rt_minutes', float),
        'rtinseconds': ('rt_seconds', float),
        'polarity': ('polarity', lambda p: p.lower()[:3]),
        'ionmode': ('polarity', lambda p: p.lower()[:3]),
        'ion': ('ion', do_nothing),
        'adduct': ('ion', do_nothing),
        'title': ('title', do_nothing),
        'ccs': ('ccs', float),
        'spectype': ('spectrum_type', do_nothing),
        'collisionenergy': ('collision_energy', do_nothing),
        'scans': ('num_scans', int),
    }

    kwargs: dict[str, int | float | str] = {}

    for line in inpt:
        # e.g. CCS=350.07224 or PEPMASS=1821.95234
        if '=' not in line:  # peak list
            continue
        key, val = line.split('=', 1)
        # strip whitespaces
        key = key.strip().lower()
        val = val.strip().strip('\n')
        label, func = keyword_to_label_and_conversion.get(key, (key, do_nothing))
        kwargs[label] = func(val)
    return kwargs


class MsSpec(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "ms_spec"

    id: Mapped[int] = mapped_column(primary_key=True)

    mz: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    ms_level: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    charge: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    rt_seconds: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    ion: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    rt_minutes: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    peaks_id: Mapped[Optional[int]] = mapped_column(ForeignKey("peak_list.id"))
    peaks: Mapped[Optional["PeakList"]] = relationship()

    feature_mgf_id: Mapped[int] = mapped_column(ForeignKey("mgf_features.id"), nullable=False)
    feature_mgf: Mapped["FeatureMgf"] = relationship(back_populates="ms_specs")


class PolarityEnum(PyEnum):
    POS = "pos"
    NEG = "neg"
    POSITIVE = "positive"
    NEGATIVE = "negative"


class FeatureMgf(SqlBaseClass, FeatureBaseClass):
    """Multiple entries can belong to a single feature id
    (when deconvolution fails or when multiple adducts are found).
    This object sets ambiguous parameters by choosing from the preferred ion
    order"""
    __tablename__ = "mgf_features"

    id: Mapped[int] = mapped_column(primary_key=True)
    feature_id: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    polarity: Mapped[Optional[str]] = mapped_column(
        Enum("pos", "neg", name="polarity_enum"),
        nullable=True
    )
    has_multiple_adducts: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)

    mz: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    charge: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    rt_seconds: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    ion: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    rt_minutes: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    ms1_id: Mapped[Optional[int]] = mapped_column(ForeignKey("peak_list.id"))
    ms1: Mapped[Optional["PeakList"]] = relationship(foreign_keys=[ms1_id])

    ms2_id: Mapped[Optional[int]] = mapped_column(ForeignKey("peak_list.id"))
    ms2: Mapped[Optional["PeakList"]] = relationship(foreign_keys=[ms2_id])

    ms_specs: Mapped[list["MsSpec"]] = relationship(
        back_populates="feature_mgf", cascade="all, delete-orphan"
    )

    combined_feature_id: Mapped[Optional[int]] = mapped_column(ForeignKey('features.id'))
    combined_feature: Mapped[Optional["FeatureCombined"]] = relationship(back_populates='mgf')

    @classmethod
    def from_lines(cls, inpt: list[str]) -> Self:
        # filter out properties that are accepted by object
        entries = {k: v for k, v in parse_ion_props(inpt) if k in cls.__annotations__}
        return cls(**entries)

    def __post_init__(self) -> None:
        """Set properties from preferred ion"""
        if (self.ms_specs is None) or len(self.ms_specs) == 0:
            self.has_multiple_adducts = False
            return

        adducts: dict[tuple[str, int], "MsSpec"] = {(msspec.ion, msspec.ms_level): msspec for msspec in self.ms_specs}
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
    # f = test()

    pl = PeakList(mzs=[1, 2, 3], intensities=[3, 4, 5])
    spec = MsSpec(mz=100, ms_level=1, charge=1, peaks=pl)

    mgf_feat = FeatureMgf(polarity='pos', ms_specs=[spec])
