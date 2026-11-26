from functools import cached_property
from typing import Optional

import pandas as pd
from sqlalchemy import Integer, Float, String, ForeignKey, Boolean
from sqlalchemy.orm import Mapped, mapped_column, relationship

from msIO.environmental.sample import Sample
from msIO.features.base import FeatureBaseClass, SqlBaseClass

METABOSCAPE_CSV_RENAME_COLUMNS: dict[str, str] = {
    'FEATURE_ID': 'feature_id',
    'RT': 'rt_seconds',
    'PEPMASS': 'M_metaboscape',
    'CCS': 'CCS',
    'SIGMA_SCORE': 'sigma_score',
    'NAME_METABOSCAPE': 'name_metaboscape',
    'MOLECULAR_FORMULA': 'formula_metaboscape',
    'ADDUCT': 'adduct_metaboscape',
    'KEGG': 'KEGG',
    'CAS': 'CAS',
    'm/z meas.': 'mz_meas',
    'Ions': 'ion_types',
    'Flags': 'flags',
    'AQ': 'annotation_type',
    'Annotation Source': 'annotation_source',
    'Include': 'include_flag'
}


def mz_from_mass(M: float | None, adduct: str | None) -> float | None:
    """
    m/z of a feature from its neutral mass and adduct (e.g. "ION=[M+H]+" as in the MetaboScape export), None if one of
    them is missing or the adduct can't be parsed.
    """
    from LipidCalculator import Adduct
    if M is None or M != M or not isinstance(adduct, str) or not adduct.strip():  # M != M: NaN
        return None
    try:
        return float(Adduct(adduct.strip()).mass_to_mz(float(M)))
    except Exception:
        return None


class Intensity(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "intensities"

    id: Mapped[int] = mapped_column(primary_key=True)

    value: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    feature_id: Mapped[int] = mapped_column(ForeignKey("metaboscape_features.feature_id"), nullable=False)
    feature: Mapped["FeatureMetaboScape"] = relationship(back_populates="intensities")

    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"), nullable=False)
    sample: Mapped["Sample"] = relationship("Sample", back_populates="intensities")


class FeatureMetaboScape(SqlBaseClass, FeatureBaseClass):
    """Container for features living in the exported feature table from MetaboScape"""
    __tablename__ = "metaboscape_features"
    # __allow_unmapped__ = True

    # TODO: define Annotation objects

    id: Mapped[int] = mapped_column(primary_key=True)  # include only if not inherited
    feature_id: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    rt_seconds: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    M_metaboscape: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    CCS: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    sigma_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    # msms_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    name_metaboscape: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    formula_metaboscape: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    adduct_metaboscape: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    KEGG: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    CAS: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    mz_meas: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    flags: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    annotation_type: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    annotation_source: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    include_flag: Mapped[bool] = mapped_column(Boolean, nullable=True)

    intensities: Mapped[list[Intensity]] = relationship(
        back_populates="feature",
        cascade="all, delete-orphan"
    )

    combined_feature_id: Mapped[Optional[int]] = mapped_column(ForeignKey('features.id'))
    combined_feature: Mapped[Optional["FeatureCombined"]] = relationship(back_populates='metaboscape')

    @classmethod
    def from_dataframe_row(cls, ser: pd.Series, sample_name_to_sample: dict[str, Sample]):
        processed = {'intensities': []}

        for k, v in ser.items():
            if k.endswith('MaxIntensity') or k.endswith('MeanIntensity'):
                continue
            if ((k not in METABOSCAPE_CSV_RENAME_COLUMNS.keys())
                    and (k not in METABOSCAPE_CSV_RENAME_COLUMNS.values())):
                # print(f'column "{k}" could not be converted, assuming it contains intensities')
                # set nan values to 0
                if isinstance(v, str):
                    v = float(v)
                v = 0. if not (v > 0) else float(v)
                if k not in sample_name_to_sample:
                    sample_name_to_sample[k] = Sample(sample_name=k)
                processed['intensities'].append(Intensity(sample=sample_name_to_sample[k], value=v))
            else:
                if k in METABOSCAPE_CSV_RENAME_COLUMNS.values():
                    k_new = k
                else:
                    k_new = METABOSCAPE_CSV_RENAME_COLUMNS[k]
                processed[k_new] = cls._convert_type(k_new, v)
        return cls(**processed)

    @classmethod
    def from_metaboscape_api(cls, feature_table_api):
        # TODO
        ...

    @property
    def M(self):
        return self.M_metaboscape

    @property
    def formula(self):
        return self.formula_metaboscape

    @cached_property
    def mz(self):
        """MetaboScape does not export mz values, so need to calculate mz from M and adduct (:"""
        from LipidCalculator import Adduct
        return Adduct(self.adduct).mass_to_mz(self.M)

    @property
    def adduct(self):
        return self.adduct_metaboscape


if __name__ == '__main__':
    pass