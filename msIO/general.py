from typing import Optional

from sqlalchemy import Integer, CheckConstraint, String
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


class Author:
    # TODO: author
    ...


class Literature:
    # TODO: literature
    ...