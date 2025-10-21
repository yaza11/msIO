"""also define back_populates"""
from __future__ import annotations

from typing import List

from sqlalchemy import Column
from sqlalchemy import Table
from sqlalchemy import ForeignKey
from sqlalchemy import Integer
from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import DeclarativeBase
from sqlalchemy.orm import relationship


class Base(DeclarativeBase):
    pass


association_table = Table(
    "association_table",
    Base.metadata,
    Column("left_id", ForeignKey("left_table.id"), primary_key=True),
    Column("right_id", ForeignKey("right_table.id"), primary_key=True),
)


class Parent(Base):
    __tablename__ = "left_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    children: Mapped[List[Child]] = relationship(
        secondary=association_table, back_populates="parents"
    )


class Child(Base):
    __tablename__ = "right_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    parents: Mapped[List[Parent]] = relationship(
        secondary=association_table, back_populates="children"
    )