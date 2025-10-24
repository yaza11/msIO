from msIO.features.base import SqlBaseClass


class Annotation(SqlBaseClass):
    """This object is a container for theoretical compounds to be used for feature and/or peak annotations.

    The structure is inspired by the form layout in the 'TargetCompounds' package.
    """
    name: str


class TargetCompound:
    """
    A target compound in the database is uniquely identified by its name. It can reference multiple peaks across
    measurements."""

    name: str
    other_names: list[str]

    samples: list['Sample']

