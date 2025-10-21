"""
(parent can have many children)

need to choose delete behavior for children

"""
from typing import List

from sqlalchemy import ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

class Base:
    pass


class Parent(Base):
    __tablename__ = "parent_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    children: Mapped[List["Child"]] = relationship()
    # or if children should be deleted if parent is deleted
    children: Mapped[List["Child"]] = relationship(
        cascade="all, delete-orphan",
        passive_deletes=True
    )


class Child(Base):
    __tablename__ = "child_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    parent_id: Mapped[int] = mapped_column(ForeignKey("parent_table.id"))

