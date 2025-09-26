from typing import Optional

from sqlalchemy import String, Float, Integer, CheckConstraint, ForeignKey, Boolean
from sqlalchemy.orm import Mapped, mapped_column, relationship

from msIO.features.base import SqlBaseClass


class Time(SqlBaseClass):
    __tablename__ = "times"

    id: Mapped[int] = mapped_column(primary_key=True)

    second: Mapped[Optional[int]] = mapped_column(
        Integer,
        CheckConstraint("second BETWEEN 0 AND 59"),
        nullable=True)
    minute: Mapped[Optional[int]] = mapped_column(
        Integer,
        CheckConstraint("minute BETWEEN 0 AND 59"),
        nullable=True)
    hour: Mapped[Optional[int]] = mapped_column(
        Integer,
        CheckConstraint("hour BETWEEN 0 AND 23"),
        nullable=True)

    time_zone: Mapped[Optional[str]] = mapped_column(
        String, nullable=True)

    day: Mapped[Optional[int]] = mapped_column(
        Integer,
        CheckConstraint("day BETWEEN 1 AND 31"),
        nullable=True)

    month: Mapped[Optional[int]] = mapped_column(
        Integer,
        CheckConstraint("month BETWEEN 1 AND 12"),
        nullable=True)

    year: Mapped[Optional[int]] = mapped_column(
        Integer,
        nullable=True
    )

    sample: Mapped[Optional["Sample"]] = relationship("Sample", back_populates='time')


class Location(SqlBaseClass):
    __tablename__ = "locations"

    id: Mapped[int] = mapped_column(primary_key=True)

    latitude: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    longitude: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    depth: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    sample: Mapped[Optional["Sample"]] = relationship("Sample", back_populates='location')


class Sample(SqlBaseClass):
    __tablename__ = "samples"

    id: Mapped[int] = mapped_column(primary_key=True)

    sample_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    # TODO: add sample code objects
    sample_code: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    sample_code_type: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    is_blank: Mapped[bool] = mapped_column(Boolean, nullable=True)
    tags: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    weight_amount: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    weight_unit: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    volume_amount: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    volume_unit: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    concentration_percent: Mapped[Optional[float]] = mapped_column(Float, CheckConstraint("concentration_percent BETWEEN 0 AND 100"), nullable=True)

    location_id: Mapped[Optional[int]] = mapped_column(ForeignKey('locations.id'))
    location: Mapped[Optional["Location"]] = relationship(back_populates='sample')

    time_id: Mapped[Optional[int]] = mapped_column(ForeignKey('times.id'))
    time: Mapped[Optional["Time"]] = relationship(back_populates='sample')

    intensities: Mapped[list["Intensity"]] = relationship("Intensity", back_populates="sample",
                                                          cascade="all, delete-orphan", passive_deletes=True)

