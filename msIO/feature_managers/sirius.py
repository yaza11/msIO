import os

import numpy as np
import pandas as pd

from msIO.feature_managers.base import FeatureManager
from msIO.features.sirius import FeatureSirius

SIRIUS_FILE_NAMES = [
    'formula_identifications',
    'compound_identifications',
    'canopus_formula_summary'
]

RENAME_FORMULA_IDENTIFICATIONS = {
    'formulaRank': 'formula_rank',
    'molecularFormula': 'formula_sirius',
    'adduct': 'adduct_sirius',
    'ZodiacScore': 'zodiac_score',
    'SiriusScore': 'sirius_score',
    'TreeScore': 'tree_score',
    'IsotopeScore': 'isotope_score',
    'numExplainedPeaks': 'num_explained_peaks',
    'explainedIntensity': 'explained_intensity',
    'medianMassErrorFragmentPeaks(ppm)': 'median_mass_error_fragments_ppm',
    'massErrorPrecursor(ppm)': 'mass_error_precursor_ppm',
    'lipidClass': 'lipid_class',
    'retentionTimeInSeconds': 'rt_seconds',
    'featureId': 'feature_id',
    'id': 'sirius_compound_folder'
}

RENAME_COMPOUND_IDENTIFICATIONS = {
    'confidenceRank': 'confidence_rank',
    'structurePerIdRank': 'structure_per_id_rank',
    'formulaRank': 'formula_rank',
    '#adducts': 'num_adducts',
    '#predictedFPs': 'num_predicted_fingerprints',
    'ConfidenceScore': 'confidence_score',
    'CSI:FingerIDScore': 'finger_id_score',
    'ZodiacScore': 'zodiac_score',
    'SiriusScore': 'sirius_score',
    'molecularFormula': 'formula_sirius',
    'adduct': 'adduct_sirius',
    'InChI': 'inchi',
    'name': 'name_sirius',
    'smiles': 'smiles',
    'xlogp': 'xlogp',
    'retentionTimeInSeconds': 'rt_seconds',
    'featureId': 'feature_id',
    'id': 'sirius_compound_folder'
}

RENAME_CANOPUS_FORMULA_SUMMARY = {
    'id': 'sirius_compound_folder',
    'molecularFormula': 'formula_sirius',
    'adduct': 'adduct_sirius',
    'NPC#pathway': 'npc_pathway_name',
    'NPC#pathway Probability': 'npc_pathway_probability',
    'NPC#superclass': 'npc_superclass_name',
    'NPC#superclass Probability': 'npc_superclass_probability',
    'NPC#class': 'npc_class_name',
    'NPC#class Probability': 'npc_class_probability',
    'ClassyFire#most specific class': 'cf_most_specific_name',
    'ClassyFire#most specific class Probability': 'cf_most_specific_probability',
    'ClassyFire#level 5': 'cf_level5_name',
    'ClassyFire#level 5 Probability': 'cf_level5_probability',
    'ClassyFire#subclass': 'cf_subclass_name',
    'ClassyFire#subclass Probability': 'cf_subclass_probability',
    'ClassyFire#class': 'cf_class_name',
    'ClassyFire#class Probability': 'cf_class_probability',
    'ClassyFire#superclass': 'cf_superclass_name',
    'ClassyFire#superclass probability': 'cf_superclass_probability',
    'ClassyFire#all classifications': 'cf_path',
    'featureId': 'feature_id',
}


def read_compound_info(path_compound_folder: str) -> dict[str, str]:
    """Read rows as key, value pairs"""
    out = {}
    path_file = os.path.join(path_compound_folder, 'compound.info')
    with open(path_file, 'r') as f:
        for line in f:
            k, v = line.split('\t')
            out[k] = v.strip('\n')
    return out


def get_sirius_file_for_tag(file: str, tag: str | None) -> str:
    return f'{file}.tsv' if tag is None else f'{file}_{tag}.tsv'


class SiriusImportManager(FeatureManager):
    _tables: dict[str, pd.DataFrame] = None
    _features: dict[int, FeatureSirius] = None

    def __init__(
            self,
            path_folder_export: str = None,
            export_tag: str = None
    ) -> None:
        self.path_folder_export = path_folder_export
        self.export_tag = export_tag
        self._read_tables()

    @property
    def feature_ids(self) -> np.ndarray:
        assert self._tables is not None, 'feature_ids not available before setting tables'
        return self._tables['formula_identifications'].feature_id.unique()

    def _read_tables(self):
        def process_with_rename(table_name: str, renamer: dict[str, str]) -> None:
            table = self._tables[table_name]
            table.rename(columns=renamer, inplace=True)
            table = table.loc[:, list(renamer.values())]
            self._tables[table_name] = table

        files = [get_sirius_file_for_tag(f, self.export_tag)
                 for f in SIRIUS_FILE_NAMES]

        self._tables: dict[str, pd.DataFrame] = {
            name: pd.read_csv(
                os.path.join(self.path_folder_export, file), sep='\t'
            )
            for name, file in zip(SIRIUS_FILE_NAMES, files)
        }

        renamers = [RENAME_FORMULA_IDENTIFICATIONS,
                    RENAME_COMPOUND_IDENTIFICATIONS,
                    RENAME_CANOPUS_FORMULA_SUMMARY]
        for n, renamer in zip(SIRIUS_FILE_NAMES, renamers):
            process_with_rename(n, renamer)

    def _inner_missing_feature(self, f_id) -> None:
        f = FeatureSirius.from_tables(f_id, self._tables)
        self._features[f_id] = f


if __name__ == '__main__':
    path_test_folder = r'\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\SIRIUS\test'
    path_full_folder = r'\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\SIRIUS\5.8.1'

    test_folder_1 = os.path.join(path_test_folder, '1_timsTOF_combined_re.sirius_1')

    props = read_compound_info(test_folder_1)

    sr = SiriusImportManager(path_folder_export=path_full_folder, export_tag='all')

    f = sr.get_feature(2)
