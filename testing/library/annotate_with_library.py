"""
This module allows to either: find annotations for a specific compound or: find features matching a specific compound.
"""
from LipidCalculator import Adduct

from msIO.feature_managers.db import FeatureManagerDB, Library
from msIO.features.metaboscape import FeatureMetaboScape

# get databases for measured data
path_file_meas = r"C:\Users\Yannick Zander\Nextcloud2\Promotion\msIO\testing\height.sqlite"
# get databases for library used to annotate the data
path_file_lib = r"C:\Users\Yannick Zander\Nextcloud2\Promotion\msIO\testing\test data\lib\archlipids_high_conf.sqlite"

# open databases
meas = FeatureManagerDB(path_file_meas)
lib = Library(path_file_lib)

# we need to specify in which uncertainty window we want to look for matches, either in mDa or ppm
max_dmz_ppm = None
max_dmz_da = 3e-3

# we need another field to specify which metric we want to use to compare MS2 similarities
metric = 'modified_cosine_greedy'

# select specific features from the measured data: use specific feature_ids for now ... in the future we will need an interactive plot to select those
f_id_for_annotation: int = 96

# only necessary for legacy data
if all([mz is None for mz in meas.mzs.values()]):
    Ms = meas.get_all_attributes_from(FeatureMetaboScape, 'M_metaboscape')
    adds = meas.get_all_attributes_from(FeatureMetaboScape, 'adduct_metaboscape')
    mzs = {f_id: Adduct(adds[f_id]).mass_to_mz(Ms[f_id]) for f_id in meas.feature_ids}
else:
    mzs = meas.mzs

# get corresponding mz value
mz_for_annotation: float = mzs[f_id_for_annotation]
# and MS1/MS2 spectrum
ms1_for_annotation = meas.get_ms_spectrum(f_id_for_annotation, level=1)
ms2_for_annotation = meas.get_ms_spectrum(f_id_for_annotation, level=2)

matched_f_ids_lib: list[int] = lib.find_matches(
    mzs={f_id_for_annotation: mz_for_annotation},
    max_dmz_ppm=max_dmz_ppm,
    max_dmz_da=max_dmz_da,
    metric=metric,
    ms2_spectra={f_id_for_annotation: ms2_for_annotation},
    return_nhits_ms2=True,
    min_ms2_score=0.3,
)

# display matched results
lib.plot_match(matched_f_ids_lib[f_id_for_annotation][0], meas, f_id_for_annotation)


# matched_f_ids_lib: list[int] = lib.find_matches(
#     mzs=mzs,
#     max_dmz_ppm=max_dmz_ppm,
#     max_dmz_da=max_dmz_da,
#     metric=metric,
#     ms2_spectra=meas.get_ms_spectra(meas.feature_ids, level=2),
#     return_nhits_ms2=True,
#     min_ms2_score=0.3,
# )