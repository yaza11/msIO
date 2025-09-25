from msIO.sql.session import get_sessionmaker
from sqlalchemy.orm import joinedload, load_only
from sqlalchemy import select, inspect
from sqlalchemy.orm import selectinload, joinedload, Load

# need to import so that sqlalchemy knows about relationships
from msIO.features.gnps import FeatureGnpsNode
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.mgf import FeatureMgf
from msIO.features.sirius import FeatureSirius
from msIO.features.base import SqlBaseClass
from msIO.features.combined import FeatureCombined


def eager_options_for(cls, strategy="selectin", maxdepth=10, seen=None):
    """ Build loader options to eagerly load all relationships on cls up to
    maxdepth.
    strategy: "selectin" (good default) or "joined".
    """
    if seen is None:
        seen = set()

    opts = []
    mapper = inspect(cls)
    for rel in mapper.relationships:
        # Skip if already visited this path (prevents cycles via backrefs)
        key = (mapper.class_, rel.key)
        if (key in seen) or key[1].startswith('combined'):
            continue
        seen_add = seen | {key}
        # Choose loader
        loader_fn = selectinload if strategy == "selectin" else joinedload
        loader = loader_fn(getattr(mapper.class_, rel.key))

        # Recurse into target mapper
        if maxdepth > 1:
            subopts = eager_options_for(rel.mapper.class_, strategy=strategy, maxdepth=maxdepth - 1, seen=seen_add)
            if subopts:
                loader = loader.options(*subopts)

        opts.append(loader)
    return opts


class FeatureManagerDB:
    """
    Object for accessing and modifying a DB created with a ProjectImportManager instance
    """

    def __init__(self, path_file_db: str):
        self.path_file = path_file_db

    @property
    def session_maker(self) -> 'session_maker':
        _Session = get_sessionmaker(self.path_file)
        return _Session

    def get_feature(self, feature_id: int) -> FeatureCombined:
        feature_id = int(feature_id)

        opts = eager_options_for(FeatureCombined, strategy="selectin")

        with self.session_maker() as session:
            obj = session.execute(
                select(FeatureCombined)
                .where(FeatureCombined.id == feature_id)
                .options(*opts)
            ).unique().scalar_one()
        return obj

    def _get_all_attributes_from(self, obj) -> list:
        with self.session_maker() as session:
            attrs = session.execute(select(obj)).scalars().all()
        return attrs

    @property
    def feature_ids(self) -> list[int]:
        return self._get_all_attributes_from(FeatureCombined.feature_id)

    @property
    def mzs(self) -> dict[int, float]:
        # cannot use properties directly
        stmt = (
            select(FeatureMetaboScape)
            .options(load_only(
                FeatureMetaboScape.feature_id,
                FeatureMetaboScape.M_metaboscape,
                FeatureMetaboScape.adduct_metaboscape
            ))
        )

        with self.session_maker() as session:
            objs = session.execute(stmt).scalars().all()
        return {o.feature_id: o.mz for o in objs}

    def find_objects_for_attr(self, attr_name: str) -> list[object]:
        raise NotImplementedError()


if __name__ == '__main__':
    from msIO.features.combined import FeatureCombined

    dbm = FeatureManagerDB(r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\SQL\database.db")

    # Session = dbm.get_active_session()
    # with Session() as session:
    #     f = session.query(FeatureCombined).filter(FeatureCombined.feature_id == 1).first()
    #     print(f.metaboscape.CCS)

    f = dbm.get_feature(1)
    print(f.metaboscape.M, f.metaboscape.adduct, f.metaboscape.mz)

    f_ids = dbm.feature_ids
    mzs = dbm.mzs
