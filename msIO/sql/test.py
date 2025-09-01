from sqlalchemy import select, or_

from msIO import PeakList
from msIO.features.base import SqlBaseClass
from msIO.features.mgf import MsSpec, FeatureMgf
from msIO.list_of_ions.base import Peak
from msIO.sql.session import get_sessionmaker, get_engine

# import objects to ensure they are in the session
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.gnps import FeatureGnpsNode
from msIO.features.sirius import FeatureSirius

from typing import Dict, List, Type


def get_by_feature_id(session, fid: int) -> Dict[Type, List]:
    """
    Return all ORM objects from ANY mapped class that have a 'feature_id'
    column equal to the given fid, grouped by class.
    """
    results: Dict[Type, List] = {}

    # Get the registry (holds all mappers) from the declarative base
    registry = SqlBaseClass.registry

    for mapper in registry.mappers:
        cls = mapper.class_
        if 'feature_id' in mapper.columns:
            matches = session.query(cls).filter(cls.feature_id == fid).all()
            if matches:
                results[cls] = matches

    return results


db_file = 'database.db'

engine = get_engine(db_file)
Session = get_sessionmaker(db_file)

# let sqlalchemy inspect database
# metadata = SqlBaseClass.metadata
# metadata.reflect(engine)


# %% filter by feature id
with Session() as session:
    feature_id_to_find = 2
    objs_by_class = get_by_feature_id(session, feature_id_to_find)

# %% filter by fragment
fragment_target_mass = 18.010565

ppm_tolerance = 5
tol = fragment_target_mass * ppm_tolerance / 1e6

stmt = (
    select(FeatureMgf.feature_id)
    .join(FeatureMgf.ms_specs)
    .join(PeakList, MsSpec.peaks_id == PeakList.id)
    .join(PeakList.peaks)
    .filter(Peak.mz.between(fragment_target_mass - tol,
                            fragment_target_mass + tol))
    .distinct()
)

with Session() as session:
    res = session.scalars(stmt).all()

