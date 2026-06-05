# create a network in which compounds with their mzs within a certain
#  similarity are connected
#   edge weights are MS2 cosine similarity scores
import numpy as np
import pandas as pd
from sqlalchemy import select
from sqlalchemy.orm import load_only

from msIO.feature_managers.db import FeatureManagerDB
from msIO.features.mgf import FeatureMgf
from msIO.features.sirius import CompoundCandidate
from msIO.sql.session import get_sessionmaker

db_file = r'\\hlabstorage.dmz.marum.de\scratch\Yannick\compounds\sql\library.sql'

dbm = FeatureManagerDB(
    db_file
)

stmt = (
    select(FeatureMgf)
    .options(load_only(
        FeatureMgf.combined_feature_id,
        FeatureMgf.mz,
    ))
)

Session = get_sessionmaker(db_file)
with Session() as session:
    objs = session.execute(stmt).scalars().all()
mzs: dict[int, float] = {o.combined_feature_id: o.mz for o in objs}

# FeatureSirius.feature_id is the same as FeatureCombined.feature_id
stmt = (
    select(CompoundCandidate)
    .options(load_only(
        CompoundCandidate.feature_id,
        CompoundCandidate.name_sirius,
    ))
)

with Session() as session:
    objs = session.execute(stmt).scalars().all()
names: dict[int, float] = {o.feature_id: o.name_sirius for o in objs}

mzs_a = np.array(list(mzs.values()))
dmzs = mzs_a[:, None] - mzs_a[None, :]

f_ids = list(mzs.keys())
dmzs_df = pd.DataFrame(data=dmzs, index=f_ids, columns=f_ids)

dmz_tol = 5e-3

mask_tri = dmzs_df.index.to_numpy()[None, :] > dmzs_df.columns.to_numpy()[:, None]
mask_dup = dmzs_df.abs().to_numpy() < dmz_tol

mask = mask_dup & mask_tri

# idcs = [(i, j) for i in range(len(f_ids)) for j in range(len(f_ids)) if mask[i, j]]

idcs = []
for i in range(len(f_ids)):
    for j in range(len(f_ids)):
        if mask[i, j]:
            idcs.append((i, j))
