import os
from functools import cached_property

import networkx as nx
import numpy as np
import pandas as pd

from msIO.feature_managers.base import FeatureManager
from msIO.features.gnps import FeatureGnpsNode, GNPS_RENAME


class GnpsImportManager(FeatureManager):
    def __init__(self, path_gnps_folder=None, path_file_gnps_graphml=None):
        assert (path_gnps_folder is not None) or (path_file_gnps_graphml is not None)
        if path_file_gnps_graphml is None:
            path_file_gnps_graphml = self._find_gnps_file(path_gnps_folder)
        self.path_file_gnps_graphml = path_file_gnps_graphml

        G: nx.Graph = nx.read_graphml(self.path_file_gnps_graphml)
        keep_columns_to_dtype = {'name': int, 'componentindex': int, 'RTConsensus': float, 'precursor mass': float}
        nodes_data = dict(map(lambda x: (int(x[0]), x[1]), G.nodes.data()))
        self._df_nodes = (
            pd.DataFrame(nodes_data)
            .loc[list(keep_columns_to_dtype), :]
            .T
            .astype(keep_columns_to_dtype)
            .rename(columns=GNPS_RENAME | dict(name='feature_id'))
            .set_index('feature_id')
        )

    @staticmethod
    def _find_gnps_file(path_gnps_folder):
        for file in os.listdir(folder := os.path.join(path_gnps_folder, 'gnps_molecular_network_graphml')):
            if file.endswith('.graphml'):
                return os.path.join(folder, file)
        raise FileNotFoundError(f'could not find graphml file in {folder}')

    def _inner_missing_feature(self, f_id) -> None:
        f = FeatureGnpsNode(
            feature_id=f_id,
            **self._df_nodes.loc[f_id, :].to_dict(),
        )
        self._features[f_id] = f

    @cached_property
    def feature_ids(self) -> np.ndarray[int]:
        return self._df_nodes.index.values


if __name__ == '__main__':
    path_gnps_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\GNPS\gnps.graphml"

    gnps = GnpsImportManager(path_file_gnps_graphml=path_gnps_file)

    f = gnps.get_feature(1)
