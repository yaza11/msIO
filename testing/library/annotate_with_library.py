"""
This module allows to either: find annotations for a specific compound or: find features matching a specific compound.
"""
import os.path

import pandas as pd
from LipidCalculator import Adduct
from sqlalchemy import select
from tqdm import tqdm

from msIO.feature_managers.db import FeatureManagerDB, Library, eager_options_for
from msIO.features.combined import FeatureCombined
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.sirius import CompoundCandidate
from msIO.sql.session import get_sessionmaker


def add_annotation_to_measured_db(f_id_meas: int, f_id_lib: int, results: pd.DataFrame, lib: Library, meas: FeatureManagerDB) -> None:
    """Add a matched library entry to the measured database."""
    ms2_score: float = results.loc[
        (results.feature_id_meas == f_id_meas) & (results.feature_id_lib == f_id_lib), 'ms2_score'
    ].squeeze()
    meas.add_annotations_from_library(lib, f_id_meas, f_id_lib, ms2_score)


# get databases for measured data
path_file_meas = r"C:\Users\Yannick Zander\Nextcloud2\Avin\database_gdgts.sqlite"
# get databases for library used to annotate the data
path_file_lib = r"C:\Users\Yannick Zander\Nextcloud2\Avin\archlipids_high_conf.sqlite"

# open databases
meas = FeatureManagerDB(path_file_meas)
lib = Library(path_file_lib)

# we need to specify in which uncertainty window we want to look for matches, either in mDa OR ppm
max_dmz_ppm = None
max_dmz_da = 3e-3

# we need another field to specify which metric we want to use to compare MS2 similarities
metric = 'modified_cosine_greedy'

# only necessary for legacy data
if all([mz is None for mz in meas.mzs.values()]):
    Ms = meas.get_all_attributes_from(FeatureMetaboScape, 'M_metaboscape')
    adds = meas.get_all_attributes_from(FeatureMetaboScape, 'adduct_metaboscape')
    mzs = {f_id: Adduct(adds[f_id]).mass_to_mz(Ms[f_id]) for f_id in meas.feature_ids}
else:
    mzs = meas.mzs

matched_f_ids_lib: dict[int, list[dict]] = lib.find_matches(
    mzs=mzs,
    max_dmz_ppm=max_dmz_ppm,
    max_dmz_da=max_dmz_da,
    metric=metric,
    ms2_spectra=meas.get_ms_spectra(meas.feature_ids, level=2),
    return_nhits_ms2=True,
    min_ms2_score=0.3,
)

# turn results into dataframe
series = []
for f_id_meas, matches in matched_f_ids_lib.items():
    for match in matches:
        series.append(pd.Series(name=f_id_meas, data=match))

results = (
    pd.concat(series, axis=1).T
    .reset_index(drop=False, names='feature_id_meas')
    .rename(columns={'feature_id': 'feature_id_lib'})
    .astype({'feature_id_meas': int, 'feature_id_lib': int, 'ms2_score': float, 'name': str, 'formula': str})
    .sort_values(by=['feature_id_meas', 'ms2_score'])
)


# select which feature from the measurement and data to plot
# (in the GUI this should be selected from the table)
f_id_meas = 6780
f_id_lib = 27192

# display matched results
lib.plot_match(f_id_lib, meas, f_id_meas)

# add annotation to measured database:
add_annotation_to_measured_db(f_id_meas, f_id_lib, results, lib, meas)