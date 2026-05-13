from typing import Any

from tqdm import tqdm

from msIO import PeakList
from msIO.sql.session import get_sessionmaker
from sqlalchemy.orm import joinedload, load_only
from sqlalchemy import select, inspect
from sqlalchemy.orm import selectinload, joinedload, Load

# need to import so that sqlalchemy knows about relationships
from msIO.features.gnps import FeatureGnpsNode
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.mgf import FeatureMgf, MsSpec
from msIO.features.sirius import FeatureSirius, CompoundCandidate, FormulaCandidate
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

    def _get_dict_for_attributes(self, parent_obj, attr) -> dict:
        vals = self._get_all_attributes_from(getattr(parent_obj, attr))
        ids = self._get_all_attributes_from(getattr(parent_obj, 'feature_id'))
        return dict(zip(ids, vals))

    def get_all_attributes_from(
            self,
            parent_obj,
            attr_name: str,
            missing_value=None,
            add_missing_values: bool = False
    ) -> dict[int, Any]:
        d = self._get_dict_for_attributes(parent_obj, attr_name)
        if add_missing_values:
            missing_features = set(self.feature_ids) - set(d.keys())
            for f in missing_features:
                d[f] = missing_value
        return d

    @property
    def molecular_mass(self) -> dict[int, float]:
        return self._get_dict_for_attributes(FeatureMetaboScape, 'M_metaboscape')

    @property
    def retention_times_in_seconds(self):
        return self._get_dict_for_attributes(FeatureMetaboScape, 'rt_seconds')

    @property
    def collisional_cross_sections(self):
        return self._get_dict_for_attributes(FeatureMetaboScape, 'CCS')

    @property
    def formula_metaboscape(self):
        return self._get_dict_for_attributes(FeatureMetaboScape, 'formula_metaboscape')

    @property
    def formula_sirius(self):
        stmt = (
            select(FormulaCandidate).filter(FormulaCandidate.formula_rank == 1)
        )

        with self.session_maker() as session:
            objs = session.execute(stmt).scalars().all()
            formulas = {}
            for o in objs:
                formulas[o.feature_id] = o.formula_sirius
        return formulas

    @property
    def name_sirius(self):
        """
        Returns dict mapping feature id to a sirius name based on the best formula. If there is no name for the
        highest scoring formula, None will be assigned to that feature."""
        stmt = (
            select(CompoundCandidate).filter(CompoundCandidate.formula_rank == 1)
        )
        names = {}
        with self.session_maker() as session:
            objs = session.execute(stmt).scalars().all()
            for o in objs:
                names[o.feature_id] = o.name_sirius
        return names

    def get_ms_spectrum(self, feature_id: int, level: int) -> PeakList:
        # fetch ms spectrum sql file for specified feature id and level (1 for isotope pattern, 2 for fragment spectrum)
        if level not in (1, 2):
            raise ValueError("level must be 1 or 2")

        stmt = (
            select(PeakList)
            .select_from(MsSpec)
            .join(MsSpec.feature_mgf)
            .join(MsSpec.peaks)
            .options(selectinload(PeakList.peaks))
            .where(
                FeatureMgf.feature_id == feature_id,
                MsSpec.ms_level == level,
                MsSpec.peaks_id.is_not(None),
            )
        )

        with self.session_maker() as session:
            return session.scalars(stmt).first()

    def find_objects_for_attr(self, attr_name: str) -> list[object]:
        raise NotImplementedError()


if __name__ == '__main__':
    from msIO.features.combined import FeatureCombined

    dbm = FeatureManagerDB(r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas really new method\SQL\database.db")

    spec = dbm.get_ms_spectrum(feature_id=8, level=2)
    spec.plot()

    # Session = dbm.get_active_session()
    # with Session() as session:
    #     f = session.query(FeatureCombined).filter(FeatureCombined.feature_id == 1).first()
    #     print(f.metaboscape.CCS)

    # f = dbm.get_feature(1)
    # print(f.metaboscape.M, f.metaboscape.adduct, f.metaboscape.mz)
    #
    # f_ids = dbm.feature_ids
    # mzs = dbm.mzs
    # Ms = dbm.molecular_mass
    # RTs = dbm.retention_times_in_seconds
    # CCS = dbm.collisional_cross_sections
    # F_m = dbm.formula_metaboscape
    F_s = dbm.formula_sirius

    # names = dbm.name_sirius

    # t = dbm.get_all_attributes_from(FeatureMetaboScape, 'M_metaboscape')
