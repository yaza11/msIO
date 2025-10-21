"""
(parent can have many children)

need to choose delete behavior for children

for bidirectional, the child to parent relationship is many to one
"""
from typing import List

from sqlalchemy import ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

class Base:
    pass


class Parent(Base):
    # this is "one"
    __tablename__ = "parent_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    children: Mapped[List["Child"]] = relationship(back_populates="parent")  # could use Set here instead of List


class Child(Base):
    # this is "many"
    __tablename__ = "child_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    parent_id: Mapped[int] = mapped_column(ForeignKey("parent_table.id"))
    parent: Mapped["Parent"] = relationship(back_populates="children")

