"""Keep this separated from main database for now"""
from typing import Optional, List

from sqlalchemy import Integer, String, ForeignKey, Boolean, Float, create_engine, Table, Column
from sqlalchemy.orm import Mapped, relationship, mapped_column, Session
from sqlalchemy.orm import DeclarativeBase

from msIO.core import PeakBaseClass


class SqlBaseClassComp(DeclarativeBase):
    pass


compound_group_to_compound_association_table = Table(
    "association_table",
    SqlBaseClassComp.metadata,
    Column("compound_group_id", ForeignKey("compound_group.id")),
    Column("compound_id", ForeignKey("compound.id")),
)


class CompoundGroup(SqlBaseClassComp):
    __tablename__ = "compound_group"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String, nullable=False)
    abbreviation: Mapped[str | None] = mapped_column(String, nullable=True)

    molecule_id: Mapped[Optional[int]] = mapped_column(ForeignKey("molecule.id"),
                                                       nullable=True)
    molecule: Mapped[Optional["Molecule"]] = relationship()

    compounds: Mapped[List['Compound']] = relationship(
        secondary=compound_group_to_compound_association_table,
        back_populates="compound_groups"
    )

    # --- Self-referential FK for parent group ---
    parent_id: Mapped[Optional[int]] = mapped_column(
        ForeignKey("compound_group.id"), nullable=True)

    # --- Relationships for hierarchy ---
    parent: Mapped[Optional["CompoundGroup"]] = relationship(
        remote_side="CompoundGroup.id",  # resolves circular relationship
        back_populates="children"
    )

    children: Mapped[List["CompoundGroup"]] = relationship(
        back_populates="parent",
        cascade="all, delete-orphan"
    )


class Compound(SqlBaseClassComp):
    __tablename__ = "compound"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String, unique=True, nullable=False)

    # many to one
    molecule_id: Mapped[Optional[int]] = mapped_column(ForeignKey("molecule.id"),
                                                       nullable=True)
    molecule: Mapped[Optional["Molecule"]] = relationship()

    # Relationship: Compound → IonPeak
    ions: Mapped[List["IonPeak"]] = relationship(
        back_populates="compound", cascade="all, delete-orphan"
    )

    compound_groups: Mapped[List["CompoundGroup"]] = relationship(
        secondary=compound_group_to_compound_association_table,
        back_populates="compounds"
    )


class Molecule(SqlBaseClassComp):
    __tablename__ = "molecule"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)

    smiles: Mapped[str | None] = mapped_column(String, nullable=True)
    inchi: Mapped[str | None] = mapped_column(String, nullable=True)
    inchkey: Mapped[str | None] = mapped_column(String, nullable=True)
    formula: Mapped[str | None] = mapped_column(String, nullable=True)
    M: Mapped[float | None] = mapped_column(Float, nullable=True)

    charge: Mapped[int | None] = mapped_column(Float, nullable=True)


class IonPeak(SqlBaseClassComp, PeakBaseClass):
    __tablename__ = "ion_peak"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    compound_id: Mapped[int] = mapped_column(ForeignKey("compound.id"), nullable=False)

    adduct_type: Mapped[Optional[str]] = mapped_column(String)
    is_ion: Mapped[bool] = mapped_column(Boolean, default=True)

    compound: Mapped["Compound"] = relationship(back_populates="ions")

    # Relationship: IonPeak → IsotopePeak
    isotopes: Mapped[List["IsotopePeak"]] = relationship(
        back_populates="ion_peak", cascade="all, delete-orphan"
    )

    # many to one
    molecule_id: Mapped[Optional[int]] = mapped_column(ForeignKey("molecule.id"),
                                                       nullable=True)
    molecule: Mapped[Optional["Molecule"]] = relationship()


class IsotopePeak(SqlBaseClassComp, PeakBaseClass):
    __tablename__ = "isotope_peak"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    ion_peak_id: Mapped[int] = mapped_column(ForeignKey("ion_peak.id"), nullable=False)

    isotope_order: Mapped[Optional[int]] = mapped_column(Integer)
    is_isotope: Mapped[bool] = mapped_column(Boolean, default=True)

    ion_peak: Mapped["IonPeak"] = relationship(back_populates="isotopes")

    # Relationship: IsotopePeak → FragmentPeak
    fragments: Mapped[List["FragmentPeak"]] = relationship(
        back_populates="isotope_peak", cascade="all, delete-orphan"
    )

    # many to one
    molecule_id: Mapped[Optional[int]] = mapped_column(ForeignKey("molecule.id"),
                                                       nullable=True)
    molecule: Mapped[Optional["Molecule"]] = relationship()


class FragmentPeak(SqlBaseClassComp, PeakBaseClass):
    __tablename__ = "fragment_peak"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)

    is_fragment: Mapped[bool | None] = mapped_column(Boolean, default=True)
    is_neutral_loss: Mapped[bool | None] = mapped_column(Boolean, nullable=True)

    isotope_peak_id: Mapped[int | None] = mapped_column(
        ForeignKey("isotope_peak.id"), nullable=False)
    isotope_peak: Mapped["IsotopePeak"] = relationship(back_populates="fragments")

    # many to one
    molecule_id: Mapped[Optional[int]] = mapped_column(ForeignKey("molecule.id"),
                                                       nullable=True)
    molecule: Mapped[Optional["Molecule"]] = relationship()


def test():
    # --- SQLITE in-memory engine ---
    engine = create_engine("sqlite:///:memory:", echo=False)

    # create tables
    SqlBaseClassComp.metadata.create_all(engine)

    # --- TEST INSERTION ---
    with Session(engine) as session:
        # Create a Compound
        ipl = CompoundGroup(name='intact polar lipids', abbreviation='IPL')
        head_1g = CompoundGroup(name='head', abbreviation='1G', parent=ipl)

        mol = Molecule(
            formula="C10H20O",
        )

        cmp = Compound(
            name="TestCompound",
            molecule=mol
        )

        cmp.compound_groups.append(ipl)
        cmp.compound_groups.append(head_1g)

        # Create IonPeak
        ion = IonPeak(
            mz=250.13,
            rt=5.5,
            intensity=1000.0,
            adduct="[M+H]+",
            adduct_type="protonated",
            is_ion=True
        )

        # Link Ion to Compound
        cmp.ions.append(ion)

        # Create IsotopePeak
        iso = IsotopePeak(
            mz=251.13,
            rt=5.5,
            intensity=500.0,
            isotope_order=0,
            is_isotope=True
        )

        ion.isotopes.append(iso)

        # Create FragmentPeak
        frag = FragmentPeak(
            mz=150.07,
            rt=5.6,
            intensity=200.0,
            is_fragment=True
        )

        iso.fragments.append(frag)

        # Commit to DB
        session.add(cmp)
        session.commit()

    # --- TEST QUERY ---
    with Session(engine) as session:
        compounds = session.query(Compound).all()
        for c in compounds:
            print(f"Compound: {c.name}, molecule={c.molecule}")
            for ion in c.ions:
                print(f"  IonPeak: mz={ion.mz}, adduct={ion.adduct}")
                for iso in ion.isotopes:
                    print(f"    IsotopePeak: mz={iso.mz}, order={iso.isotope_order}")
                    for frag in iso.fragments:
                        print(f"      FragmentPeak: mz={frag.mz}, intensity={frag.intensity}")


if __name__ == '__main__':
    test()
