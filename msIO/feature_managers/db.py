from typing import Any, Iterable, Literal, Callable

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from msIO import PeakList
from msIO.environmental.sample import Sample
from msIO.metrics import cosine_similarity_sim
from msIO.sql.session import get_sessionmaker
from sqlalchemy.orm import load_only
from sqlalchemy import select, inspect
from sqlalchemy.orm import selectinload, joinedload

# need to import so that sqlalchemy knows about relationships
from msIO.features.metaboscape import FeatureMetaboScape, Intensity
from msIO.features.mgf import FeatureMgf, MsSpec
from msIO.features.sirius import CompoundCandidate, FormulaCandidate
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
    Object for accessing and modifying a DB created with a ProjectImportManager
     instance
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
                # FeatureMetaboScape.M_metaboscape,
                # FeatureMetaboScape.adduct_metaboscape
                FeatureMetaboScape.mz_meas,
            ))
        )

        with self.session_maker() as session:
            objs = session.execute(stmt).scalars().all()
        return {o.feature_id: o.mz_meas for o in objs}

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
        return self._get_dict_for_attributes(
            FeatureMetaboScape, 'formula_metaboscape'
        )

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

    @property
    def names_metaboscape(self) -> dict[int, str]:
        return self.get_all_attributes_from(FeatureMetaboScape, 'name_metaboscape')

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

    def _get_ms_spectra_limited_variable_number(
            self,
            feature_ids: list[int],
            level: int
    ) -> dict[int, PeakList]:
        if len(feature_ids) == 0:
            return {}

        stmt = (
            select(FeatureMgf.feature_id, PeakList)
            .select_from(MsSpec)
            .join(FeatureMgf, MsSpec.feature_mgf_id == FeatureMgf.id)
            .join(PeakList, MsSpec.peaks_id == PeakList.id)
            .options(selectinload(PeakList.peaks))
            .where(
                FeatureMgf.feature_id.in_(feature_ids),
                MsSpec.ms_level == level,
                MsSpec.peaks_id.is_not(None),
            )
        )

        with self.session_maker() as session:
            rows = session.execute(stmt).all()

        return {
            feature_id: peak_list
            for feature_id, peak_list in rows
        }

    def get_ms_spectra(
            self,
            feature_ids: list[int],
            level: int,
            max_spectra_per_query: int = 20_000
    ) -> dict[int, PeakList]:
        level = int(level)

        if level not in (1, 2):
            raise ValueError("level must be 1 or 2")

        # split feature_ids into chunks of max_spectra_per_query
        feature_id_chunks: list[list[int]] = [
            feature_ids[i:i + max_spectra_per_query]
            for i in range(0, len(feature_ids), max_spectra_per_query)
        ]

        out: dict[int, PeakList] = {}
        for feature_id_chunk in feature_id_chunks:
            print('num f_ids:', len(feature_id_chunk))
            out |= self._get_ms_spectra_limited_variable_number(
                feature_id_chunk, level
            )
            print('num spectra after:', len(out))

        return out

    def get_intensities(self, feature_id) -> dict[str, int]:
        """
        Returns a dict mapping sample names to intensities for the given
        feature id.
        """
        stmt = (
            select(Sample.sample_name, Intensity.value)
            .join(Intensity, Intensity.sample_id == Sample.id)
            .where(Intensity.feature_id == feature_id)
        )

        with self.session_maker() as session:
            ints = session.execute(stmt).all()
            return dict(ints)

    def compare_features(self, f_id1: int, f_id2: int):
        fig, axs = plt.subplots(nrows=4)

        f_1 = self.get_feature(f_id1)
        f_2 = self.get_feature(f_id2)

        ms1: dict[int, PeakList] = {ms_spec.ms_level: ms_spec.peaks for ms_spec in f_1.mgf.ms_specs}
        ms2: dict[int, PeakList] = {ms_spec.ms_level: ms_spec.peaks for ms_spec in f_2.mgf.ms_specs}

        # MS/MS mirror plot
        if 2 in ms1:
            axs[0].stem(ms1[2].mzs, ms1[2].intensities, label=f_id1, markerfmt='', linefmt='C0-')
        if 2 in ms2:
            axs[0].stem(ms2[2].mzs, -np.asarray(ms2[2].intensities), label=f_id2, markerfmt='', linefmt='C1-')
        axs[0].set_title('MS/MS spectra')
        axs[0].legend()
        axs[0].set_xlabel('m/z in Da')

        # MS1 plot
        if 1 in ms1:
            axs[1].stem(np.asarray(ms1[1].mzs) - ms1[1].mzs[0],
                        np.asarray(ms1[1].intensities) / max(ms1[1].intensities),
                        label=f_id1, markerfmt='', linefmt='C0-')
        if 1 in ms2:
            axs[1].stem(np.asarray(ms2[1].mzs) - ms2[1].mzs[0], np.asarray(ms2[1].intensities) / max(ms2[1].intensities),
                        label=f_id2, markerfmt='', linefmt='C1--')
        axs[1].set_title('To M0 shifted and 1-scaled MS1 spectra')
        axs[1].legend()
        axs[1].set_xlabel('m/z in Da')

        # table with properties
        df = pd.DataFrame(dict(
            rt_min=[round(f_1.metaboscape.rt_seconds / 60, 2), round(f_2.metaboscape.rt_seconds / 60, 2)],
            ccs=[round(f_1.metaboscape.CCS, 1), round(f_2.metaboscape.CCS, 1)],
            mz=[round(f_1.metaboscape.mz_meas, 4), round(f_2.metaboscape.mz_meas, 4)],
            main_ion=[f_1.metaboscape.adduct, f_2.metaboscape.adduct],
            M=[round(f_1.metaboscape.M_metaboscape, 4), round(f_1.metaboscape.M_metaboscape, 4)],
            annotation=[f_1.metaboscape.name_metaboscape, f_2.metaboscape.name_metaboscape],
        ))
        df.index = [f_id1, f_id2]
        axs[2].axis("off")
        pd.plotting.table(ax=axs[2], data=df.T, loc="center", cellLoc="center", edges='open')

        # bar plot with intensities
        ints1 = {i.sample.sample_name: i.value for i in f_1.metaboscape.intensities}
        ints2 = {i.sample.sample_name: i.value for i in f_2.metaboscape.intensities}
        df = pd.DataFrame(data=[ints1, ints2], index=[f_id1, f_id2]).T

        df.plot.bar(rot=45, ax=axs[3])
        axs[3].set_title('Intensities')
        axs[3].legend()

        return fig, axs

    def find_objects_for_attr(self, attr_name: str) -> list[object]:
        raise NotImplementedError()


class Library(FeatureManagerDB):
    """
    This class is not intended to be structured like this longterm. This is
    merely necessary because the current architecture is used for turning msp
    libraries into sql files.
    """
    _mzs: dict[int, float] = None
    _names: dict[int, str] = None

    _mzs_sorted: np.ndarray[float] = None
    _f_ids_sorted: np.ndarray[int] = None

    @property
    def names(self) -> dict[int, str]:
        stmt = (
            select(CompoundCandidate)
        )
        names = {}
        with self.session_maker() as session:
            objs = session.execute(stmt).scalars().all()
            for o in objs:
                names[o.feature_id] = o.name_sirius
        return names

    def _set_sorted_mzs(self):
        _mzs: np.ndarray[float] = np.asarray(list(self.mzs.values()))
        o = np.argsort(_mzs)
        self._mzs_sorted = _mzs[o]
        self._f_ids_sorted: np.ndarray[int] = np.asarray(list(self.mzs.keys()))[o]

    @property
    def f_ids_sorted(self) -> np.ndarray[int]:
        if self._f_ids_sorted is None:
            self._set_sorted_mzs()
        return self._f_ids_sorted

    @property
    def mzs_sorted(self) -> np.ndarray[float]:
        if self._mzs_sorted is None:
            self._set_sorted_mzs()
        return self._mzs_sorted

    def find_matches_precursor(
            self,
            mzs: Iterable[float],
            max_dmz_da: float = None,
            max_dmz_ppm: float | int = None,
    ) -> list[list[int]]:
        """Returns the matched feature ids"""
        assert (max_dmz_da is None) ^ (max_dmz_ppm is None), \
            'provide either max_dmz_da or max_dmz_ppm (but not both)'

        mzs = np.asarray(mzs)

        if max_dmz_da is None:
            max_dmz_da = mzs * (max_dmz_ppm * 1e-6)

        # get all features that match the mz within the tolerance
        idcs_left = np.searchsorted(self.mzs_sorted, mzs - max_dmz_da, side='right')
        idcs_right = np.searchsorted(self.mzs_sorted, mzs + max_dmz_da, side='right')

        f_ids: list[np.ndarray[int]] = [
            self.f_ids_sorted[idx_left:idx_right]
            for idx_left, idx_right in zip(idcs_left, idcs_right)
        ]

        # convert types
        return [[int(f_id) for f_id in _f_ids] for _f_ids in f_ids]

    def find_matches(
            self,
            mzs: Iterable[float],
            max_dmz_da: float=None,
            max_dmz_ppm: float | int = None,
            ms2: Iterable[PeakList] = None,
            max_ms2_dmz_da: float = 10e-3,
            min_sim: float = 0.7,
            metric: Callable[[PeakList | None, PeakList | None], float] | Literal['cosine_fwd', 'cosine_bwd', 'cosine_sim'] = 'cosine_sim',
    ):
        assert (max_dmz_da is None) ^ (max_dmz_ppm is None), \
            'provide either max_dmz_da or max_dmz_ppm (but not both)'

        f_ids_matched_precursor: list[np.ndarray[int]] = self.find_matches_precursor(
            mzs, max_dmz_da, max_dmz_ppm
        )

        # pre-load ms2 spectra of all matched features
        f_ids_matched_unique: set[int] = set()
        for _f_ids_matched_precursor in f_ids_matched_precursor:
            f_ids_matched_unique.update(_f_ids_matched_precursor)

        matches_ms2: dict[int, PeakList] = self.get_ms_spectra(
            list(f_ids_matched_unique), level=2
        )

        # score
        ...

        # return features above thr
        ...

        # # retrieve the feature ids and load the complete objects
        # options = []
        #
        # # Only load MS2 spectra + peaks if actually needed
        # if ms2 is not None:
        #     options.append(
        #         selectinload(
        #             FeatureMgf.ms_specs.and_(MsSpec.ms_level == 2)
        #         )
        #         .selectinload(MsSpec.peaks)
        #         .selectinload(PeakList.peaks)
        #     )
        #
        # stmt = (
        #     select(FeatureMgf, CompoundCandidate)
        #     .outerjoin(
        #         CompoundCandidate,
        #         CompoundCandidate.feature_id == FeatureMgf.combined_feature_id
        #     )
        #     .where(FeatureMgf.mz.between(mz_lower, mz_upper))
        # )
        #
        # with self.session_maker() as session:
        #     matches = session.execute(stmt).unique().all()

        def add_notnone_attributes(obj):
            return {
                k: v
                for k, v in obj.__dict__.items()
                if (v is not None) and (not k.startswith('_'))
            }

        matches_transformed: list[dict] = []
        for match in matches:
            feature_mgf, compound_candidate = match
            # all features that are not None
            out = add_notnone_attributes(compound_candidate) | add_notnone_attributes(feature_mgf)

            if ms2 is not None:
                ms2_refs: list[PeakList] = [
                    ms_spec.peaks
                    for ms_spec in feature_mgf.ms_specs
                    if ms_spec.ms_level == 2
                ]
                if len(ms2_refs) == 0:  # no ms2 information
                    continue
                score = cosine_similarity_sim(ms2_refs[0], ms2, max_ms2_dmz_da)
                if score < min_cos_sim:
                    continue
            else:
                score = None

            matches_transformed.append(out)
            out['dmz'] = feature_mgf.mz - mz
            if score is not None:
                out['cosine_score'] = score

        return matches_transformed


if __name__ == '__main__':
    # lib_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\compounds\sql\library.sql"
    # lib_file = r"C:\Users\yanni\Downloads\library_complete.sql"
    # meas_file = r"C:\Users\yanni\Downloads\Guaymas new method height recursive\SQL\database.db"
    lib_file = r"C:\Users\Yannick Zander\Downloads\library_complete.sql"
    meas_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\SQL\database.db"


    # target_mzs = [653.68090, 667.69624, 802.74909, 1073.80953] * 100
    # target_names = ['Archaeol(20:0_20:0)', 'MeO-Archaeol(20:0_20:0)', 'Rib-Archaeol(20:0_20:0)', 'SQ-Archaeol(25:3_30:5)'] * 100

    lib = Library(lib_file)
    mzs_lib = lib.mzs
    names_lib = lib.names

    meas = FeatureManagerDB(meas_file)
    target_mzs = list(meas.mzs.values())
    names_metaboscape = meas.names_metaboscape
    target_names = [names_metaboscape[f_id] for f_id in meas.mzs]

    matched_f_ids: list[list[int]] = lib.find_matches_precursor(target_mzs, max_dmz_da=5e-3)
    matched_names = [[names_lib[f_id] for f_id in f_ids] for f_ids in matched_f_ids]

    self = lib

    f_ids_matched_precursor: list[list[int]] = self.find_matches_precursor(target_mzs, max_dmz_da=5e-3)

    # pre-load ms2 spectra of all matched features
    f_ids_matched_unique: set[int] = set()
    for _f_ids_matched_precursor in f_ids_matched_precursor:
        f_ids_matched_unique.update(_f_ids_matched_precursor)

    matches_ms2: dict[int, PeakList] = self.get_ms_spectra(list(f_ids_matched_unique), level=2)

    meas_f_with_match = [f_id for f_id, matched_f_id in zip(meas.feature_ids, matched_f_ids) if len(matched_f_id) > 0]
    meas_ms2: dict[int, PeakList] = meas.get_ms_spectra(meas_f_with_match, level=2)

    # assign scores
    ms2_scores = [
        [cosine_similarity_sim(matches_ms2.get(f_id_lib), meas_ms2.get(f_id_meas), 10e-3) for f_id_lib in f_id_libs]
        for f_id_meas, f_id_libs in zip(meas.feature_ids, matched_f_ids)
    ]

    import time
    t0 = time.time()

    # ms2_dummy = PeakList(mzs=[100, 200, 300], intensities=[0.1, 0.2, 0.3])
    #
    # for target_mz, target_name in tqdm(zip(target_mzs, target_names), total=len(target_mzs)):
    #     matches = lib.find_matches(target_mz, max_dmz_da=5e-3, min_cos_sim=0)
    #     for match in matches:
    #         if match['name_sirius'].strip() == target_name:
    #             break
    #     else:
    #         print(f'target name {target_name} is not in the matches')

    t1 = time.time()

    print(f'{t1 - t0:.2f} seconds')

    # match = matches[0]



    # from msIO.features.combined import FeatureCombined
    #
    # dbm = FeatureManagerDB(r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\SQL\database.db")
    #
    # ints = dbm.get_intensities(feature_id=1)
    #
    #
    # self = dbm
    # f_id1, f_id2 = 6000, 6050
    # self.compare_features(f_id1, f_id2)
    #
    # spec = dbm.get_ms_spectrum(feature_id=6050, level=2)
    # spec.plot()
    #
    # specs = dbm.get_ms_spectra(feature_ids=list(range(1, 1000)), level=2)

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
    # F_s = dbm.formula_sirius

    # names = dbm.name_sirius

    # t = dbm.get_all_attributes_from(FeatureMetaboScape, 'M_metaboscape')
