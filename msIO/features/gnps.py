import json

from dataclasses import dataclass
from typing import Self, Optional

from sqlalchemy import Integer, Float, String, ForeignKey
from sqlalchemy.orm import mapped_column, Mapped, relationship

from msIO.features.base import FeatureBaseClass, SqlBaseClass

GNPS_RENAME: dict[str, str] = {
    'componentindex': 'cluster_label',
    'precursor mass': 'M_gnps',
    'RTConsensus': 'rt_seconds'
}


class FeatureGnpsNode(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "gnps_features"

    id: Mapped[int] = mapped_column(primary_key=True)  # required unless inherited
    feature_id: Mapped[Optional[int]] = mapped_column(Integer)
    cluster_label: Mapped[Optional[int]] = mapped_column(Integer)
    M_gnps: Mapped[Optional[float]] = mapped_column(Float)
    rt_seconds: Mapped[Optional[float]] = mapped_column(Float)

    # not mapped to a column
    # TODO: convert to sql object as well
    #  converting to json for now
    other: Mapped[Optional[str]] = mapped_column(String)

    combined_feature_id: Mapped[int] = mapped_column(ForeignKey('features.id'))
    combined_feature: Mapped["FeatureCombined"] = relationship(back_populates='gnps')

    @classmethod
    def from_graphml(cls, inpt: tuple[str | int, dict]) -> Self:
        """Handle input from G.nodes.data()"""
        _id, props = inpt
        processed = {'feature_id': int(_id)}
        _other = {}
        for k, v in props.items():
            if k in GNPS_RENAME:
                k_renamed = GNPS_RENAME[k]
                processed[k_renamed] = cls._convert_type(k_renamed, v)
            else:
                _other[k] = v
        processed['other'] = json.dumps(_other)
        return cls(**processed)

    @property
    def M(self):
        return self.M_gnps

    @property
    def other_dict(self):
        return json.loads(self.other)


if __name__ == '__main__':
    import networkx as nx

    path_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\GNPS\gnps_molecular_network_graphml\0315b5f68fc74173b4f18a6906a0c1f9.graphml"

    G: nx.Graph = nx.read_graphml(path_file)
    fs = list(G.nodes.data())
    f = FeatureGnpsNode.from_graphml(fs[0])
