from typing import Optional

from sqlalchemy import String, Float, CheckConstraint, ForeignKey, Boolean
from sqlalchemy.orm import Mapped, mapped_column, relationship

from msIO.features.base import SqlBaseClass


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

    concentration_percent: Mapped[Optional[float]] = mapped_column(
        Float,
        CheckConstraint("concentration_percent BETWEEN 0 AND 100"),
        nullable=True
    )

    # many to one
    location_id: Mapped[Optional[int]] = mapped_column(ForeignKey('locations.id'))
    location: Mapped[Optional["Location"]] = relationship(back_populates='sample')

    # one to one
    time_id: Mapped[Optional[int]] = mapped_column(ForeignKey('times.id'))
    time: Mapped[Optional["Time"]] = relationship(back_populates='sample')

    intensities: Mapped[list["Intensity"]] = relationship("Intensity",
                                                          back_populates="sample",
                                                          cascade="all, delete-orphan",
                                                          passive_deletes=True)

    # possible belongs to multiple measurements
    measurements: Mapped[list["Measurement"]] = relationship('Measurement',
                                                             back_populates="samples",
                                                             cascade="all, delete-orphan",
                                                             passive_deletes=True)
