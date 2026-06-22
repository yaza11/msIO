from functools import cached_property
from typing import Any, Iterable, Literal, Callable

import numpy as np
import pandas as pd
import logging

from LipidCalculator.rdkit.plotting import mplt_mol
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Contrib.Glare.glare import Library
from tqdm import tqdm

from msIO import PeakList
from msIO.environmental.sample import Sample
from msIO.list_of_ions.read_mca import MoleculeAnnotation
from msIO.metrics import cosine_similarity_sym, cosine_similarity_forward, cosine_similarity_backward
from msIO.sql.session import get_sessionmaker
from sqlalchemy.orm import load_only
from sqlalchemy import select, inspect
from sqlalchemy.orm import selectinload, joinedload

# need to import so that sqlalchemy knows about relationships
from msIO.features.metaboscape import FeatureMetaboScape, Intensity
from msIO.features.mgf import FeatureMgf, MsSpec
from msIO.features.sirius import CompoundCandidate, FormulaCandidate
from msIO.features.combined import FeatureCombined

logger = logging.getLogger(__name__)


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

    @cached_property
    def feature_ids(self) -> list[int]:
        return self._get_all_attributes_from(FeatureCombined.feature_id)

    @cached_property
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

    @cached_property
    def molecular_mass(self) -> dict[int, float]:
        return self._get_dict_for_attributes(FeatureMetaboScape, 'M_metaboscape')

    @cached_property
    def retention_times_in_seconds(self):
        return self._get_dict_for_attributes(FeatureMetaboScape, 'rt_seconds')

    @cached_property
    def collisional_cross_sections(self):
        return self._get_dict_for_attributes(FeatureMetaboScape, 'CCS')

    @cached_property
    def formula_metaboscape(self):
        return self._get_dict_for_attributes(
            FeatureMetaboScape, 'formula_metaboscape'
        )

    @cached_property
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

    @cached_property
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

    @cached_property
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
            .join(MsSpec.feature_mgf)
            .join(MsSpec.peaks)
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
            out |= self._get_ms_spectra_limited_variable_number(
                feature_id_chunk, level
            )

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

    def add_annotations_from_library(self, library: "Library" = None, library_file: str | None = None):
        assert (library is None) ^ (library_file is None), 'provide either library or library_file'
        if library is not None:
            library = Library(library_file)
        ...


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

    @cached_property
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

    @cached_property
    def smiles(self) -> dict[int, str]:
        return self.get_all_attributes_from(CompoundCandidate, 'smiles')

    @cached_property
    def inchis(self) -> dict[int, str]:
        return self.get_all_attributes_from(CompoundCandidate, 'inchi')

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
            mzs: Iterable[float] | dict[int, float],
            max_dmz_da: float = None,
            max_dmz_ppm: float | int = None,
            ms2_spectra: Iterable[PeakList | None] | dict[int, PeakList | None] = None,
            max_ms2_dmz_da: float = 10e-3,
            min_ms2_score: float | None = 0.7,
            metric: Callable[[PeakList | None, PeakList | None], float] | Literal['cosine_fwd', 'cosine_bwd', 'cosine_sim'] = 'cosine_sim',
            return_nhits_ms2: bool = False,
            require_ms2: bool = False,
    ):
        assert (max_dmz_da is None) ^ (max_dmz_ppm is None), \
            'provide either max_dmz_da or max_dmz_ppm (but not both)'
        if (as_dicts := isinstance(mzs, dict)) and (ms2_spectra is not None) and (not isinstance(ms2_spectra, dict)):
            raise ValueError('If mzs are provided as dict, ms2 must also be a dict')

        if as_dicts:
            mz_ids: list[int] = list(mzs.keys())
            mzs: list[float] = list(mzs.values())
            if ms2_spectra is not None:
                ms2_spectra: list[PeakList | None] = [ms2_spectra.get(f_id) for f_id in mz_ids]
        else:
            assert len(ms2_spectra) == len(mzs), \
                'ms2 and mzs must have the same length (can set ms2 to None for some mzs, if not available)'
            mz_ids = None

        logger.info(f'finding precursor matches')
        matched_f_ids: list[list[int]] = self.find_matches_precursor(
            mzs=mzs, max_dmz_da=max_dmz_da, max_dmz_ppm=max_dmz_ppm
        )
        logger.info(f'found matches for {sum(1 for mchs in matched_f_ids if len(mchs) > 0):_} features')

        # no ms2 spectra provided, so we can not score them
        if ms2_spectra is None:
            if as_dicts:
                out = {mz_id: matches_per_meas for mz_id, matches_per_meas in zip(mz_ids, matched_f_ids) if
                       len(matches_per_meas) > 0}
                return out
            else:
                return matched_f_ids

        if isinstance(metric, str):
            if metric == 'cosine_sim':
                metric = cosine_similarity_sym
            elif metric == 'cosine_fwd':
                metric = cosine_similarity_forward
            elif metric == 'cosine_backward':
                metric = cosine_similarity_backward
            else:
                raise ValueError(f'Unknown metric {metric}')

        # fetch ms2 spectra of library matches
        matched_f_ids_raveled: set[int] = set()
        for _f_ids_matched_precursor in matched_f_ids:
            matched_f_ids_raveled.update(_f_ids_matched_precursor)

        logger.info(f'loading lib ms2 spectra for {len(matched_f_ids_raveled):_} features')
        matched_ms2_spectra_lib: dict[int, PeakList] = self.get_ms_spectra(list(matched_f_ids_raveled), level=2)

        ann_libs: dict[int, str] = self._get_dict_for_attributes(FeatureMetaboScape, 'annotation_type')

        # assign scores
        # ms2_scores: list[list[float]] = []
        # matched_f_ids_filtered: list[list[int]] = []
        matches: list[list[dict]] = []
        for mz_meas, ms2_measured, f_id_libs in tqdm(
                zip(mzs, ms2_spectra, matched_f_ids),
                desc='assigning ms2 scores',
                total=len(mzs)
        ):
            matches_per_meas: list[dict] = []
            for f_id_lib in f_id_libs:
                ms2_score, *n_hits_ms2 = metric(matched_ms2_spectra_lib.get(f_id_lib), ms2_measured, max_ms2_dmz_da, return_nhits=return_nhits_ms2)
                # only add the entry if ms2 is not required or score is above the threshold
                if (require_ms2 and np.isnan(ms2_score)) or (ms2_score < min_ms2_score):
                    continue

                mz_lib = self.mzs[f_id_lib]
                match: dict = dict(
                    feature_id=f_id_lib,
                    name=self.names.get(f_id_lib),
                    formula=self.formula_metaboscape.get(f_id_lib),
                    ms2_score=ms2_score,
                    dmz_mda= (dmz := (mz_meas - mz_lib)) * 1e3,
                    dmz_ppm= dmz / mz_lib * 1e6,
                    source_library=ann_libs.get(f_id_lib, 'unknown'),
                )
                if return_nhits_ms2:
                    match['n_hits_ms2'] = n_hits_ms2[0]
                matches_per_meas.append(match)

            matches.append(matches_per_meas)

        if as_dicts:
            out = {mz_id: matches_per_meas for mz_id, matches_per_meas in zip(mz_ids, matches) if len(matches_per_meas) > 0}
            return out
        return matches

    def plot_compound_overview(self, f_id, axs: tuple[plt.Axes, plt.Axes] = None, **kwargs):
        if axs is None:
            _, ax = plt.subplots(nrows=2)

        if self.smiles.get(f_id) is not None:
            mol = Chem.MolFromSmiles(self.smiles[f_id])
        elif self.inchis.get(f_id) is not None:
            mol = Chem.MolFromInchi(self.inchis[f_id])
        else:
            mol = Chem.Mol()
        mplt_mol(mol=mol, ax=axs[0], **kwargs)
        axs[0].set_title(f'feature id: {f_id}, name: {self.names.get(f_id, "")}, mz: {self.mzs[f_id]:.4f} Da')
        ms2: PeakList = self.get_ms_spectrum(f_id, level=2)
        ms2.plot(ax=axs[1])
        return axs

    def compare_compounds(self, f_ids):
        fig, axs = plt.subplots(nrows=2, ncols=len(f_ids))

        for i, f_id in enumerate(f_ids):
            self.plot_compound_overview(f_id, axs=axs[:, i])
        return fig, axs



if __name__ == '__main__':
    import time
    # lib_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\compounds\sql\library.sql"
    # lib_file = r"C:\Users\yanni\Downloads\library_complete.sql"
    # meas_file = r"C:\Users\yanni\Downloads\Guaymas new method height recursive\SQL\database.db"
    lib_file = r"C:\Users\Yannick Zander\Downloads\library_complete.sqlite"
    # lib_file = r"C:\Users\Yannick Zander\Downloads\library_ipl.sqlite"

    lib = Library(lib_file)
    lib.compare_compounds([334520, 481732, 686932, 457245, 334430, 334431])
    # lib.plot_compound_overview(427291)
    # lib.plot_compound_overview(427275)
    # lib.plot_compound_overview(427276)

    # meas_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\SQL\database.db"
    #
    # # target_mzs = [653.68090, 667.69624, 802.74909, 1073.80953] * 100
    # # target_names = ['Archaeol(20:0_20:0)', 'MeO-Archaeol(20:0_20:0)', 'Rib-Archaeol(20:0_20:0)', 'SQ-Archaeol(25:3_30:5)'] * 100
    #
    # logging.basicConfig(level=logging.INFO)
    #
    # t0 = time.time()
    # print('loading library')
    # lib = Library(lib_file)
    #
    # print('loading measurement')
    # meas = FeatureManagerDB(meas_file)
    # mzs_meas: dict[int, float] = meas.mzs
    # ms2_meas = meas.get_ms_spectra(meas.feature_ids, level=2)
    #
    # # matches = find_matches(
    # #     lib,
    # #     mzs_meas,
    # #     ms2_spectra=ms2_meas,
    # #     max_ms2_dmz_da=10e-3,
    # #     min_ms2_score=.7,
    # #     max_dmz_da=5e-3,
    # #     require_ms2=False
    # # )
    # matches: dict[int, list[dict]] = lib.find_matches(
    #     mzs_meas,
    #     ms2_spectra=ms2_meas,
    #     max_ms2_dmz_da=10e-3,
    #     min_ms2_score=.7,
    #     max_dmz_da=5e-3,
    #     require_ms2=False,
    #     return_nhits_ms2=True
    # )
    # t1 = time.time()
    #
    # print(f'finding annotations took {(t1 - t0) / 60:.2f} minutes')
    #
    # def plot_matches(f_id_meas):
    #     matched_anns: list[dict] = matches[f_id_meas]
    #
    #     fig, axs = plt.subplots(nrows=len(matched_anns) + 1, ncols=1, sharex=True)
    #
    #     ax=axs[0]
    #     ms2_meas[f_id_meas].plot(ax)
    #     ax.set_title(f'name MetaboScape: {meas.names_metaboscape[f_id_meas]}')
    #
    #     for i, matched_ann in enumerate(matched_anns):
    #         ax = axs[i+1]
    #         ms = lib.get_ms_spectrum(matched_ann['feature_id'], level=2)
    #         if ms is None:
    #             continue
    #         ms.plot(ax)
    #         ax.set_title(f'name: {matched_ann['name']}, ppm: {matched_ann["dmz_ppm"]:.1f}, MS score: {matched_ann["ms2_score"]:.2f}')
    #     return fig, axs
    # #
    # # def to_mca() -> list[MoleculeAnnotation]:
    # #     ...
    # #
    # plt.close()
    # plot_matches(4719)
    # plt.show()


    # assign certainty levels
    ...


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
