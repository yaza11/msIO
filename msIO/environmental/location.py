from typing import Optional

from sqlalchemy import String, Float
from sqlalchemy.orm import Mapped, mapped_column, relationship

from msIO.environmental.sample import Sample
from msIO.features.base import SqlBaseClass


class Location(SqlBaseClass):
    __tablename__ = "locations"

    id: Mapped[int] = mapped_column(primary_key=True)

    name: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    latitude: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    longitude: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    depth: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    sample: Mapped[Optional["Sample"]] = relationship("Sample", back_populates='location')
