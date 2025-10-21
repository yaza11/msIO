from typing import Optional

from sqlalchemy import ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, relationship

class Base:
    pass


class Parent(Base):
    __tablename__ = "parent_table"

    id: Mapped[int] = mapped_column(primary_key=True)
    # place foreign key in parent table
    child_id: Mapped[Optional[int]] = mapped_column(ForeignKey("child_table.id"))
    child: Mapped[Optional["Child"]] = relationship()


class Child(Base):
    __tablename__ = "child_table"

    id: Mapped[int] = mapped_column(primary_key=True)