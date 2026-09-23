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

annotations_to_add: dict[int, int] = {f_id_meas: f_id_lib}

# add annotation to measured database:
Session = get_sessionmaker(path_file_meas)

with Session() as session:
    opts = eager_options_for(FeatureCombined, strategy="selectin")
    for f_id_meas, f_id_lib in tqdm(annotations_to_add.items(), desc='adding features to DB', total=len(annotations_to_add)):
        f_meas = session.execute(
                select(FeatureCombined)
                .where(FeatureCombined.feature_id == f_id_meas)
                .options(*opts)
            ).unique().scalar_one()

        match_properties: pd.Series = results.loc[
            (results.feature_id_meas == f_id_meas) & (results.feature_id_lib == f_id_lib)
            , :
        ].squeeze()

        # add attributes from library match to measured feature
        f_lib: FeatureCombined = lib.get_feature(f_id_lib)  # read only
        # add as compound candidate
        compound_candidate_lib: CompoundCandidate = f_lib.sirius.compound_candidates[0]
        compound_candidate_meas: CompoundCandidate = CompoundCandidate(
            num_adducts=1,
            confidence_score=match_properties['ms2_score'],
            adduct_sirius=f_lib.metaboscape.adduct_metaboscape,
            name_sirius=compound_candidate_lib.name_sirius,
            smiles=compound_candidate_lib.smiles,
            inchi=compound_candidate_lib.inchi,
            xlogp=compound_candidate_lib.xlogp,
            sirius_compound_folder=f'manual_library_annotation from "{os.path.basename(path_file_lib)}" for feature "{f_id_lib}"'
        )
        f_meas.sirius.compound_candidates.append(compound_candidate_meas)
    session.commit()
