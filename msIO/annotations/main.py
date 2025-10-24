from typing import Optional, List

from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import Mapped, relationship
from sqlalchemy.testing.schema import mapped_column

from msIO.features.base import SqlBaseClass


class Annotation(SqlBaseClass):
    """This object is a container for theoretical compounds to be used for feature and/or peak annotations.

    The structure is inspired by the form layout in the 'TargetCompounds' package.
    """
    name: str


class Alias(SqlBaseClass):
    __tablename__ = "aliases"

    id: Mapped[int] = mapped_column(primary_key=True)
    target_compound_id: Mapped[int] = mapped_column(ForeignKey("target_compounds.id"))

    name: Mapped[str] = mapped_column(String, nullable=False)



class Adduct:
    pass


class Compound:
    """
    A target compound in the database is uniquely identified by its name. It can reference multiple peaks across
    measurements."""

    name: str
    # one to many
    other_names: Mapped[Optional[List['Alias']]] = relationship()

    # one to many
    ms1: Spectrum
    ms2: Spectrum
    rt: float
    mz: float
    ccs: float

    adducts: List[Adduct]
    # peaks: Mapped[Optional[List['Peak']]] = relationship()

