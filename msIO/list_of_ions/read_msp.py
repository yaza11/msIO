import numpy as np
import pandas as pd
from tqdm import tqdm

from msIO.environmental.location import Location
from msIO.features.combined import FeatureCombined
from msIO.list_of_ions.base import BaseLib, PeakList

msp_key_to_py: dict[str, str] = {
    'NAME': 'name',
    'Name': 'name',
    'MW': 'mz',
    'PRECURSORMZ': 'mz',
    'PRECURSORTYPE': 'ion',
    'FORMULA': 'formula',
    'Ontology': 'ontology',
    'ONTOLOGY': 'ontology',
    'INCHIKEY': 'inchikey',
    'IONIZATION': 'ionization_method',
    'INSTRUMENTTYPE': 'instrument_type',
    'INSTRUMENTYPE': 'instrument_type',
    'INSTRUMENT': 'instrument_type',
    'SMILES': 'smiles',
    'RETENTIONTIME': 'rt_minutes',
    'CCS': 'ccs',
    'IONMODE': 'polarity',
    'COLLISIONENERGY': 'collision_energy',
    'Comment': 'comment',
    'COMMENT': 'comment',
    'Num Peaks': 'num_peaks'
}


def try_float(val):
    """attempt to convert to float, set to nan otherwise"""
    try:
        cval = float(val)
    except ValueError:
        cval = np.nan
    return cval


def float_keep_str(val):
    """keep at str if not float convertable"""
    try:
        cval = float(val)
    except ValueError:
        cval = val.strip('\n')
    return cval


def default_str(val):
    return val.strip('\n')


converts = {
    'mz': float,
    'rt_minutes': try_float,
    'ccs': try_float,
    'collision_energy': float_keep_str,
    'num_peaks': int,
    'polarity': lambda x: x.lower()
}


def _parse_lines(lines: list[str]) -> dict[str, str | int | float]:
    e = {}
    for l in lines:
        if ':' not in l:
            continue
        k, v = l.split(':', 1)
        pyk = msp_key_to_py[k]
        v = converts.get(pyk, default_str)(v)
        e[pyk] = v
    return e


class MSPReader(BaseLib):
    def __init__(self, path_lib, splitter_peaks_list='\t'):
        entries = {}
        self.peak_lists: dict[int, PeakList] = {}
        with open(path_lib, 'r') as f:
            lines = []
            for i, l in tqdm(enumerate(f)):
                lines.append(l)
                if l == '\n':
                    entries[i] = _parse_lines(lines)
                    self.peak_lists[i] = PeakList.from_lines(
                        lines, splitter=splitter_peaks_list)
                    lines = []

        self.df_features: pd.DataFrame = pd.DataFrame.from_dict(
            entries, orient='index')
        self.df_features.loc[:, 'ms_level'] = 2
        if 'rt_minutes' in self.df_features.columns:
            self.df_features.loc[:, 'rt_seconds'] = self.df_features.rt_minutes * 60

    def get_ms2(
            self,
            mz: float = None,
            rt_minutes: float = None,
            rt_seconds: float = None,
            mass_tolerance: float = 3e-3,
            rt_minutes_tolerance: float = 10 / 6,
            rt_seconds_tolerance: float = 10,
            ion: str = None,
            feature_index: int = None
    ) -> tuple[pd.DataFrame, list[PeakList]]:
        return self._get_ms2(
            mz, rt_minutes, rt_seconds, mass_tolerance,
            rt_minutes_tolerance, rt_seconds_tolerance,
            ion, feature_index
        )


def composition_msdial(msdial):
    ft = msdial.df_features
    o = {k.strip(): 0 for k in ft.ontology.unique()}
    n = set()
    for i, row in tqdm(ft.iterrows(), total=ft.shape[0]):
        if (name := row.loc['name'].strip()) not in n:
            n.add(name)
            o[row.ontology.strip()] += 1


if __name__ == '__main__':
    # path_lib = r"C:\Users\Yannick Zander\Downloads\MSMS-Public_experimentspectra-pos-VS19.msp"
    # MSDialLib = MSPReader(path_lib)
    # s = MSDialLib.get_ms2(mz=636.53323, mass_tolerance=10e-3)

    path_file = '\\\\hlabstorage.dmz.marum.de\\scratch\\Yannick\\compounds\\1G-AEG_pos.msp'

    rdr = MSPReader(path_file, splitter_peaks_list=' ')
