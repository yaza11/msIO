from LipidCalculator import CompoundDict

from msIO.feature_managers.db import FeatureManagerDB, Library, lib_match_result_to_table

path_file_lib = r"C:\Users\yanni\Nextcloud2\Avin\archlipids_high_conf.sqlite"
path_file_meas = r"C:\Users\yanni\Nextcloud2\Avin\database_gdgts_with_annotations.sqlite"


meas = FeatureManagerDB(path_file_meas)
try:
    meas.add_neutral_losses()
except ValueError as e:
    print(e)

lib = Library(path_file_lib)

mz_tol_da = 3e-3

matches = lib.find_matches(
    meas.mzs,
    max_dmz_da=mz_tol_da,
    return_nhits_ms2=True,
    metric='modified_cosine_greedy',
    ms2_spectra=meas.get_ms_spectra(meas.feature_ids, level=2),
    require_ms2=True,
    min_ms2_score=0.6
)

matches_df = lib_match_result_to_table(matches)
# add best match for each feature so that we have something more interesting to look at
matches_df.sort_values(by=['feature_id_meas', 'ms2_score'], inplace=True, ascending=[True, False])
best_matches = matches_df.groupby(by='feature_id_meas').first().reset_index()

# remove all previous annotations
anns = set(meas.get_f_ids_manual_annotations().values())

for f_id_meas, fid_lib in anns:
    meas.remove_annotation(f_id_meas, fid_lib)

for _, row in best_matches.iterrows():
    meas.add_annotations_from_library(
        library=lib,
        f_id_meas=row['feature_id_meas'],
        f_id_lib=row['feature_id_lib'],
        ms2_score=row['ms2_score']
    )

# mark all features that have a neutral loss indicating NH4 adducts
nh3_loss = CompoundDict('NH3').mass
f_ids = meas.filter_by_neutral_loss(neutral_loss=nh3_loss, mz_tol_da=.005)
meas.plot_scatter(add_annotations=True, highlight_features=f_ids)
