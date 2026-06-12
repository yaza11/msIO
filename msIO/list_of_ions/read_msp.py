import numpy as np
import pandas as pd
from tqdm import tqdm
from itertools import islice

from msIO.environmental.location import Location
from msIO.features.combined import FeatureCombined
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.mgf import FeatureMgf, MsSpec
from msIO.features.sirius import FeatureSirius, CompoundCandidate
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
        pyk = msp_key_to_py.get(k, k)
        v = converts.get(pyk, default_str)(v)
        e[pyk] = v
    return e


class MSPReader(BaseLib):
    splitter_peaks_list: str = None
    path_file: str = None
    low_memory: bool = None
    n_lines: int = None
    _current_file_offset = 0

    peak_list: dict[int, PeakList] = None

    def __init__(self, path_file=None, splitter_peaks_list=None, low_memory=False):
        self.path_file = path_file

        if splitter_peaks_list is not None:
            self.splitter_peaks_list = splitter_peaks_list

        self.low_memory = low_memory

        with open(self.path_file, 'rb') as f:
            self.n_lines = sum(1 for _ in f)
        with open(self.path_file, 'r', errors='replace', encoding='utf-8') as f:
            self.n_features = sum([1 for l in f if l == '\n']) + 1

        if path_file is None:
            return

        if not low_memory:
            self.read_file()
            # update number of features (multiple blank rows are counted as multiple features)
            self.n_features = self.df_features.shape[0]

    def _process_lines(self, lines: list[str]) -> tuple[dict, PeakList]:
        # determine splitter for peaks
        if self.splitter_peaks_list is None:
            # determine from last line
            lines_peaks = [l for l in lines if l[0].isnumeric()]
            if len(lines_peaks) == 0:
                self.splitter_peaks_list = None
            else:
                if '\t' in lines_peaks[0]:
                    self.splitter_peaks_list = '\t'
                elif ',' in lines_peaks[0]:
                    self.splitter_peaks_list = ','
                elif ' ' in lines_peaks[0]:
                    self.splitter_peaks_list = ' '
                else:
                    raise ValueError(
                        f'Unable to determine splitter from line '
                        f'{lines_peaks[0]}, please specify manually'
                    )
        entries = _parse_lines(lines)
        peak_list = PeakList.from_lines(lines, splitter=self.splitter_peaks_list)
        return entries, peak_list

    def read_next(self):
        """always assigns index 0"""
        with open(self.path_file, 'rb') as f:
            # print(self._current_file_offset)
            f.seek(self._current_file_offset)
            lines = []

            while True:
                line = f.readline()
                self._current_file_offset += len(line)
                line = line.decode('utf-8', errors="replace").strip('\r\n')
                if len(line) == 0:
                    break
                lines.append(line)

        ent, peak_list = self._process_lines(lines)
        entries = {0: ent}

        self.peak_list: dict[int, PeakList] = {0: peak_list}

        if len(ent) == 0:
            self.df_features = pd.DataFrame()
            return

        self.df_features: pd.DataFrame = pd.DataFrame.from_dict(
            entries, orient='index')
        self.df_features.loc[:, 'ms_level'] = 2
        if 'rt_minutes' in self.df_features.columns:
            self.df_features.loc[:, 'rt_seconds'] = self.df_features.rt_minutes * 60


    def read_file(self):
        entries = {}
        self.peak_list: dict[int, PeakList] = {}

        with open(self.path_file, 'r', encoding='utf-8', errors='replace') as f:
            lines = []
            for i, l in tqdm(
                    enumerate(f),
                    total=self.n_lines,
                    desc='parsing msp file',
                    smoothing=1/50
            ):
                lines.append(l)
                if not l == '\n':  # feature ends
                    continue

                ent, peak_list = self._process_lines(lines)
                entries[i] = ent
                self.peak_list[i] = peak_list

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

    def create_feature(self, idx: int, feature_id=None) -> FeatureCombined:
        """Convert row and peak list to msIO features that can be stored as sql."""
        def get_attr_or_none(attr_name):
            val = row.get(attr_name)
            if val is None:
                return None
            if (isinstance(val, int | float)) and (val < 0):
                return None
            return val

        row: dict = {k.lower(): v for k, v in self.df_features.loc[idx, :].items()}

        if idx in self.peak_list:
            peaks: PeakList = self.peak_list.get(idx)
            ms_level = row.get('ms_level')
            if ms_level is None:
                ms_level = 2
            ms_spec = MsSpec(peaks=peaks, ms_level=int(ms_level))
            ms_specs = [ms_spec]
        else:
            ms_specs = None

        # make sure we use all fields from msp to initialize objects
        f_mgf = FeatureMgf(
            feature_id=feature_id,
            ms_specs=ms_specs
        )
        f_metabo = FeatureMetaboScape(
            feature_id=feature_id,
            rt_seconds=get_attr_or_none('rt_seconds'),
            CCS = get_attr_or_none('ccs'),
            mz_meas=get_attr_or_none('mz'),
            adduct_metaboscape=get_attr_or_none('ion'),
            formula_metaboscape=get_attr_or_none('formula'),
            # abuse annotation source for comment
            annotation_source=get_attr_or_none('comment'),
        )

        compound_candidate = CompoundCandidate(
            feature_id=feature_id,
            name_sirius=get_attr_or_none('name'),
            xlogp = get_attr_or_none('logp'),
            inchi = get_attr_or_none('inchi'),
            smiles = get_attr_or_none('smiles'),
            confidence_rank = get_attr_or_none('confidence_level')  # higher rank/level is better
        )
        f_sirius = FeatureSirius(
            feature_id=feature_id,
            compound_candidates=[compound_candidate]
        )

        f_combined = FeatureCombined(
            feature_id=feature_id,
            metaboscape=f_metabo,
            mgf=f_mgf,
            sirius=f_sirius,
        )
        return f_combined


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

    # path_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\compounds\julius\fragments\1G-AEG_pos.msp"
    # path_file = r"C:\Users\yanni\Downloads\Archlips_Full_spectral_library.msp"

    path_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\compounds\julius\fragments\1G-AEG_pos.msp"

    # rdr = MSPReader(path_file, low_memory=True)
    # rdr.read_file()

    rdr = MSPReader(path_file, low_memory=False)
    # for i in tqdm(range(rdr2.n_features)):
    #     rdr2.read_next()
    # print(rdr2.df_features.T)

    f = rdr.create_feature(rdr.df_features.index[0])
