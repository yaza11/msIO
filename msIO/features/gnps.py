from dataclasses import dataclass
from typing import Self

from msIO.features.base import FeatureBaseClass

GNPS_RENAME: dict[str, str] = {
    'componentindex': 'cluster_label',
    'precursor mass': 'M_gnps',
    'RTConsensus': 'rt_seconds'
}


@dataclass
class FeatureGnpsNode(FeatureBaseClass):
    feature_id: int
    cluster_label: int
    M_gnps: float
    rt_seconds: float
    other: dict = None

    @classmethod
    def from_graphml(cls, inpt: tuple[str | int, dict]) -> Self:
        """Handle input from G.nodes.data()"""
        _id, props = inpt
        processed = {'feature_id': int(_id), 'other': {}}
        for k, v in props.items():
            if k in GNPS_RENAME:
                k_renamed = GNPS_RENAME[k]
                processed[k_renamed] = cls._convert_type(k_renamed, v)
            else:
                processed['other'][k] = v
        return cls(**processed)

    @property
    def M(self):
        return self.M_gnps


if __name__ == '__main__':
    import networkx as nx

    path_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\GNPS\gnps_molecular_network_graphml\0315b5f68fc74173b4f18a6906a0c1f9.graphml"

    G: nx.Graph = nx.read_graphml(path_file)
    fs = list(G.nodes.data())
    f = FeatureGnpsNode.from_graphml(fs[0])
