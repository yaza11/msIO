import os
from functools import cached_property
from typing import Iterable

import networkx as nx
import numpy as np
from tqdm import tqdm

from msIO.feature_managers.base import FeatureManager
from msIO.features.gnps import FeatureGnpsNode


class GnpsImportManager(FeatureManager):
    def __init__(self, path_gnps_folder=None, path_file_gnps_graphml=None):
        assert (path_gnps_folder is not None) or (path_file_gnps_graphml is not None)
        if path_file_gnps_graphml is None:
            path_file_gnps_graphml = self._find_gnps_file(path_gnps_folder)
        self.path_file_gnps_graphml = path_file_gnps_graphml

        self._nodes: dict[int, dict] = self._get_nodes()

    @staticmethod
    def _find_gnps_file(path_gnps_folder):
        for file in os.listdir(folder := os.path.join(path_gnps_folder, 'gnps_molecular_network_graphml')):
            if file.endswith('.graphml'):
                return os.path.join(folder, file)
        raise FileNotFoundError(f'could not find graphml file in {folder}')

    def _get_nodes(self) -> dict[int, dict]:
        G: nx.Graph = nx.read_graphml(self.path_file_gnps_graphml)
        return dict(map(lambda x: (int(x[0]), x[1]), G.nodes.data()))

    def _inner_missing_feature(self, f_id) -> None:
        f = FeatureGnpsNode.from_graphml((f_id, self._nodes[f_id]))
        self._features[f_id] = f

    @cached_property
    def feature_ids(self) -> np.ndarray[int]:
        return np.array(sorted(self._get_nodes().keys()))


if __name__ == '__main__':
    path_gnps_folder = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\GNPS"

    gnps = GnpsImportManager(path_gnps_folder=path_gnps_folder)

