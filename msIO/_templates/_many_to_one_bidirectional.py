"""this is the same as one to many bidirectional with the names of Parent and child swapped"""
from typing import List

from sqlalchemy import ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

class Base:
    pass


class Parent(Base):
    __tablename__ = "parent_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    child_id: Mapped[int] = mapped_column(ForeignKey("child_table.id"))
    child: Mapped["Child"] = relationship(back_populates="parents")


class Child(Base):
    __tablename__ = "child_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    parents: Mapped[List["Parent"]] = relationship(back_populates="child")