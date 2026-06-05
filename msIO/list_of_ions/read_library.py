import io

import pandas as pd
from tqdm import tqdm

from msIO import MSPReader
from msIO.list_of_ions.base import BaseLib, PeakList


rename_key = {
    'Name': 'name',
    'Formula': 'formula',
    'InChI': 'inchi',
    'InChIKey': 'inchikey',
    'Structure': 'structure',
    'CommentSpec': 'comment_spec',
    'InstType': 'instrument_type',
    'InstName': 'instrument_name',
    'IoniMethod': 'ionization_method',
    'IonPolarity': 'polarity',
    'MSMS': 'ms_level',
    'PreIon': 'mz',
    'ColEnergy': 'collision_energy',
    'Num Peaks': 'num_peaks',
    'Smiles': 'smiles'
}


class MetaboLibraryReader(BaseLib):
    def __init__(self, path_lib):
        entries = {}
        self.peak_list: dict[int, PeakList] = {}

        with open(path_lib, 'rb') as f:
            n_lines = sum(1 for _ in f)

        with open(path_lib, 'r', encoding='utf-8') as f:
            props = {}
            ints = []
            mzs = []
            for i, l in tqdm(
                    enumerate(f),
                    total=n_lines,
                    desc='parsing msp file',
                    smoothing=1 / 50
            ):
                l = l.strip('\n').rstrip(' ')
                if l == '':  # terminates entry
                    entries[i] = props
                    self.peak_list[i] = PeakList(mzs=mzs, intensities=ints)
                    props = {}
                elif l[0].isnumeric():
                    ints_and_mzs = l.split(' ')
                    for i, int_or_mz in enumerate(ints_and_mzs):
                        if (i % 2) == 0:
                            ints.append(float(int_or_mz))
                        else:
                            mzs.append(float(int_or_mz))
                else:
                    key, value = l.split(':')
                    key_renamed = rename_key.get(key, key)
                    props[key_renamed] = value

        self.df_features: pd.DataFrame = pd.DataFrame.from_dict(
            entries, orient='index')
        self.df_features.loc[:, 'ms_level'] = 2
        if 'rt_minutes' in self.df_features.columns:
            self.df_features.loc[:, 'rt_seconds'] = self.df_features.rt_minutes * 60

    def to_msp(self) -> MSPReader:
        msp = MSPReader()
        msp.peak_list = self.peak_list.copy()
        msp.df_features = self.df_features.copy()
        return msp


if __name__ == '__main__':
    path_file = r"C:\Users\yanni\Downloads\archaeol library JSL_2023_v2.library"
    reader = MetaboLibraryReader(path_file).to_msp()
