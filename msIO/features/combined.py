"""
Module for combining information from MetaboScape, SIRIUS and GNPS exports.

Each feature is associated with molecular properties, MS1, MS2, annotations, ...
"""
from typing import Literal, Iterable

from msIO import PeakList
from msIO.features.base import FeatureBaseClass
from msIO.features.gnps import FeatureGnpsNode
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.mgf import FeatureMgf
from msIO.features.sirius import FeatureSirius


class FeatureCombined(FeatureBaseClass):
    """Combine information from """

    # sirius.mgf
    feature_id: int = None
    mz: float = None
    polarity: Literal['pos', 'neg'] = None

    rt_seconds: float = None
    rt_minutes: float = None
    ms1: PeakList = None
    ms2: PeakList = None

    # metaboscape .csv
    M_metaboscape: float = None
    name_metaboscape: str = None
    formula_metaboscape: str = None
    adduct_metaboscape: str = None
    intensities: dict[str, int] = None
    CCS: float = None
    sigma_score: float = None
    CAS: float = None
    KEGG: float = None

    # gnps
    M_gnps: float
    cluster_label: int = None
    other: dict = None

    # sirius files
    #  formula
    formula_sirius: str = None
    formula_rank: int = None
    adduct_sirius: str = None
    zodiac_score: float = None
    sirius_score: float = None
    tree_score: float = None
    isotope_score: float = None
    num_explained_peaks: int = None
    explained_intensity: float = None
    median_mass_error_fragments_ppm: float = None
    mass_error_precursor_ppm: float = None
    lipid_class: str = None
    sirius_compound_folder: str = None

    #  compound candidate
    confidence_rank: int = None
    structure_per_id_rank: int = None
    num_adducts: int = None
    num_predicted_fingerprints: int = None
    confidence_score: float = None
    finger_id_score: float = None
    inchi: str = None
    name_sirius: str = None
    smiles: str = None
    xlogp: float = None

    #  compound group
    npc_pathway_name: str = None
    npc_pathway_probability: float = None
    npc_superclass_name: str = None
    npc_superclass_probability: float = None
    npc_class_name: str = None
    npc_class_probability: float = None
    cf_most_specific_name: str = None
    cf_most_specific_probability: float = None
    cf_level5_name: str = None
    cf_level5_probability: float = None
    cf_subclass_name: str = None
    cf_subclass_probability: float = None
    cf_class_name: str = None
    cf_class_probability: float = None
    cf_superclass_name: str = None
    cf_superclass_probability: float = None
    cf_path: str = None

    def __init__(self, features: Iterable[FeatureBaseClass]):
        _possible_conflicts = ['adduct', 'formula', 'rt_seconds']
        # TODO: SIRIUS should overwrite adduct and formula, rt_seconds should be from metaboscape
        features_dict = {}
        feature_types = ['metaboscape', 'mgf', 'gnps', 'sirius']
        feature_names = [FeatureMetaboScape.__name__, FeatureMgf.__name__, FeatureGnpsNode.__name__, FeatureSirius.__name__]
        feature_name_to_key = dict(zip(feature_names, feature_types))
        for feature in features:
            # find right key
            key = feature_name_to_key[feature.__class__.__name__]
            assert key not in features_dict, 'cannot have same feature type twice'
            features_dict[key] = feature

        attrs = {}
        # order matters, take less derived properties, if possible
        for f_type in ['sirius', 'gnps', 'mgf', 'metaboscape']:
            if f_type not in features_dict:
                continue
            f = features_dict[f_type]
            attrs |= {k: v for k, v in f.__dict__.items() if v is not None}

        self.__dict__ |= attrs





def load_database(in_file: str) -> dict[int, FeatureCombined]:
    ...


class FeatureSelector:
    """Quickly obtain all feature properties for a given feature ID"""
    def __init__(self):
        ...

    def __call__(self, feature_id: ...):
        ...


if __name__ == '__main__':
    pass
    # path_mgf_sirius = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.sirius.mgf"
    # path_gnps_folder = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\GNPS"
    # path_sirius_folder = r'\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\SIRIUS\5.8.1'
    #
    # mgf = MgfImportManager(path_mgf_sirius)
    # gnps = GnpsImportManager(path_gnps_folder=path_gnps_folder)
    # sr = SiriusImportManager.from_export(path_folder_export=path_sirius_folder, export_tag='all')


