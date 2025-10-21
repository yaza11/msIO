"""same as one to many bidirectional but we use Mapped[Child] instead of Mapped[List[Child]]"""
from sqlalchemy import ForeignKey, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship


class Base:
    pass

class Parent(Base):
    __tablename__ = "parent_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    child: Mapped["Child"] = relationship(back_populates="parent")


class Child(Base):
    __tablename__ = "child_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    parent_id: Mapped[int] = mapped_column(ForeignKey("parent_table.id"))
    parent: Mapped["Parent"] = relationship(back_populates="child")
    # we can enforece single parent like this (otherwise it will not be detected)
    parent: Mapped["Parent"] = relationship(back_populates="child", single_parent = True)

    # should always use this
    __table_args__ = (UniqueConstraint("parent_id"),)