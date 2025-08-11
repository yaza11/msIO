from dataclasses import dataclass
from typing import Optional

import pandas as pd
from sqlalchemy import Integer, Float, String, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

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

}


class Intensity(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "intensities"

    id: Mapped[int] = mapped_column(primary_key=True)

    sample_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    value: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)

    feature_id: Mapped[int] = mapped_column(ForeignKey("feature_metaboscape.id"), nullable=False)
    feature: Mapped["FeatureMetaboScape"] = relationship(back_populates="intensities")


@dataclass
class FeatureMetaboScape(SqlBaseClass, FeatureBaseClass):
    """Container for features living in the exported feature table from MetaboScape"""
    __tablename__ = "feature_metaboscape"
    # __allow_unmapped__ = True

    id: Mapped[int] = mapped_column(primary_key=True)  # include only if not inherited
    feature_id: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    rt_seconds: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    M_metaboscape: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    CCS: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    sigma_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    name_metaboscape: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    formula_metaboscape: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    adduct_metaboscape: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    KEGG: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    CAS: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    intensities: Mapped[list[Intensity]] = relationship(
        back_populates="feature",
        cascade="all, delete-orphan"
    )

    @classmethod
    def from_dataframe_row(cls, ser: pd.Series):
        processed = {'intensities': []}
        for k, v in ser.items():
            if k.endswith('MaxIntensity') or k.endswith('MeanIntensity'):
                continue
            if k not in METABOSCAPE_CSV_RENAME_COLUMNS.keys():
                # print(f'column "{k}" could not be converted, assuming it contains intensities')
                # set nan values to 0
                if isinstance(v, str):
                    v = float(v)
                v = 0 if not (v > 0) else int(v)
                processed['intensities'].append(Intensity(sample_name=k, value=v))
                print(f'adding intensity {v} for sample {k}')
            else:
                k_new = METABOSCAPE_CSV_RENAME_COLUMNS[k]
                processed[k_new] = cls._convert_type(k_new, v)
        print(processed['intensities'])
        return cls(**processed)

    @property
    def M(self):
        return self.M_metaboscape

    @property
    def formula(self):
        return self.formula_metaboscape

    @property
    def adduct(self):
        return self.adduct_metaboscape


if __name__ == '__main__':
    t = FeatureMetaboScape.py_types()

    path_file_csv = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.csv"

    df = pd.read_csv(path_file_csv)
    ser = df.iloc[0, :]

    f = FeatureMetaboScape.from_dataframe_row(ser)
