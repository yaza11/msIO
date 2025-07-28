from dataclasses import dataclass

import numpy as np
import pandas as pd

from msIO.features.base import FeatureBaseClass

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


@dataclass
class FeatureMetaboScape(FeatureBaseClass):
    """Container for features living in the exported feature table from MetaboScape"""
    feature_id: int
    rt_seconds: float = None
    M_metaboscape: float = None
    CCS: float = None
    sigma_score: float = None
    name_metaboscape: str = None
    formula_metaboscape: str = None
    adduct_metaboscape: str = None
    intensities: dict[str, int] = None
    KEGG: float = None
    CAS: float = None

    @classmethod
    def from_dataframe_row(cls, ser: pd.Series):
        processed = {'intensities': {}}
        for k, v in ser.items():
            if k not in METABOSCAPE_CSV_RENAME_COLUMNS.values():
                # set nan values to 0
                processed['intensities'][k] = 0 if not (v > 0) else int(v)
            else:
                processed[k] = cls._convert_type(k, v)
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
    path_file_csv = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.csv"

    df = pd.read_csv(path_file_csv)
    ser = df.iloc[0, :]

    f = FeatureMetaboScape.from_dataframe_row(ser)
