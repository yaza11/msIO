import os.path
import typing
from functools import cached_property
from typing import Iterable, Literal, Any

import numpy as np
import pandas as pd
from tqdm import tqdm

from msIO import MgfImportManager
from msIO.feature_managers.base import FeatureManager
from msIO.feature_managers.gnps import GnpsImportManager
from msIO.feature_managers.metaboscape import MetaboscapeImportManager
from msIO.feature_managers.sirius import SiriusImportManager
from msIO.features.base import SqlBaseClass
from msIO.features.combined import FeatureCombined, FeatureCombinedFlat


class ProjectImportManager(FeatureManager):
    """Object for creating a db file from import managers"""
    def __init__(
            self,
            metaboscape_manager: MetaboscapeImportManager = None,
            mgf_manager: MgfImportManager = None,
            gnps_manager: GnpsImportManager = None,
            sirius_manager: SiriusImportManager = None
    ) -> None:
        self._metaboscape_manager = metaboscape_manager
        self._mgf_manager = mgf_manager
        self._gnps_manager = gnps_manager
        self._sirius_manager = sirius_manager

        self._name_to_manager = {
            'metaboscape': self._metaboscape_manager,
            'mgf': self._mgf_manager,
            'gnps': self._gnps_manager,
            'sirius': self._sirius_manager
        }

    @property
    def active_managers(self) -> dict[str, FeatureManager]:
        return {name: manager
                for name, manager in self._name_to_manager.items()
                if manager is not None}

    @cached_property
    def feature_ids(self) -> np.ndarray[int]:
        f_ids = set()
        for m in self.active_managers.values():
            f_ids |= set(m.feature_ids)
        return np.array(sorted(f_ids))

    def _inner_missing_feature(self, f_id) -> None:
        features = {name: m.get_feature(f_id) for name, m in self.active_managers.items() if f_id in m.feature_ids}
        f = FeatureCombined(feature_id=f_id, **features)
        self._features[f_id] = f

    def get_metaboscape_feature(self, f_id: int):
        assert 'metaboscape' in self.active_managers
        return self.active_managers['metaboscape'].get_feature(f_id)

    def get_mgf_feature(self, f_id: int):
        assert 'mgf' in self.active_managers
        return self.active_managers['mgf'].get_feature(f_id)

    def get_gnps_feature(self, f_id: int):
        assert 'gnps' in self.active_managers
        return self.active_managers['gnps'].get_feature(f_id)

    def get_sirius_feature(self, f_id: int):
        assert 'sirius' in self.active_managers
        return self.active_managers['sirius'].get_feature(f_id)

    def to_sql(self, db_file: str, feature_ids=None) -> None:
        from msIO.sql.session import get_sessionmaker

        if feature_ids is None:
            feature_ids = self.feature_ids

        Session = get_sessionmaker(db_file)

        with Session() as session:
            for f_id in tqdm(feature_ids, desc='adding features to DB', total=len(feature_ids)):
                # features = [m.get_feature(f_id) for m in self.active_managers.values() if f_id in m.feature_ids]
                # for f in features:
                #     # print(f'adding {f.__class__.__name__} to session')
                #     session.add(f)
                f = self.get_feature(f_id)
                session.add(f)
            session.commit()


def write_table(project_import_manager: ProjectImportManager) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Compile feature attributes in dataframes"""
    columns = []
    excluded_columns = []

    f_ids = project_import_manager.feature_ids

    include_dtypes = [int, float, str]
    for col, dtype in FeatureCombinedFlat.__annotations__.items():
        if col in columns:
            continue
        # exclude complex dtypes
        if dtype not in include_dtypes:
            excluded_columns.append(col)
        else:
            columns.append(col)
    proto_intensities: dict[int, dict] = {}  # don't know columns for intensities beforehand
    excluded_properties: dict[int, dict[str, Any]] = {}
    df = pd.DataFrame(columns=columns, index=f_ids)
    for f_id in tqdm(f_ids, desc='writing tables', total=f_ids.shape[0]):
        f = project_import_manager.get_feature(f_id)
        for col, val in f.__dict__.items():
            if col == 'intensities':
                proto_intensities[f_id] = val
            elif col in columns:
                df.at[f_id, col] = val
            else:
                if f_id not in excluded_properties:
                    excluded_properties[f_id] = {}
                excluded_properties[f_id][col] = val

    df_intensities = pd.DataFrame(proto_intensities)

    return df, df_intensities, excluded_properties


if __name__ == '__main__':
    path_metaboscape_csv = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.csv"
    path_mgf_sirius = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape\timsTOF_combined_re.sirius.mgf"
    path_gnps_folder = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\GNPS"
    path_sirius_folder = r'\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\SIRIUS\5.8.1'

    metaboscape = MetaboscapeImportManager(path_metaboscape_csv)
    mgf = MgfImportManager(path_mgf_sirius)
    # mgf = None
    gnps = GnpsImportManager(path_gnps_folder=path_gnps_folder)
    sr = SiriusImportManager(path_folder_export=path_sirius_folder, export_tag='all')

    project_import_manager = ProjectImportManager(mgf_manager=mgf,
                                                  gnps_manager=gnps,
                                                  sirius_manager=sr,
                                                  metaboscape_manager=metaboscape)

    res = write_table(project_import_manager)

    # # %%
    # f_id = 6
    #
    # # f = project_import_manager.get_feature(f_id)
    #
    # f_mgf = project_import_manager.get_mgf_feature(f_id)
    # f_gnps = project_import_manager.get_gnps_feature(f_id)
    # f_sirius = project_import_manager.get_sirius_feature(f_id)
    # f_meta = project_import_manager.get_metaboscape_feature(f_id)
    #
    # f = FeatureCombined(feature_id=f_id, gnps=f_gnps)
    #
    # f_comb = project_import_manager.get_feature(f_id)
    # # import matplotlib.pyplot as plt
    # # f_mgf.ms2.plot(); plt.show()
    #
    # # %%
    #
    #
    # from msIO.sql.session import initiate_db
    #
    # db_file = 'database.db'
    #
    # initiate_db(db_file)
    # # project_import_manager.to_sql(db_file, feature_ids=project_import_manager.feature_ids[:10])
    # project_import_manager.to_sql(db_file)
