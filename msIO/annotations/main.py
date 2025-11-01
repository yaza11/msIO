"""Keep this separated from main database for now"""
from typing import Optional, List

from sqlalchemy import Integer, String, ForeignKey, Boolean, Float, create_engine
from sqlalchemy.orm import Mapped, relationship, mapped_column, Session
from sqlalchemy.orm import DeclarativeBase

from msIO.core import PeakBaseClass


class SqlBaseClassComp(DeclarativeBase):
    pass


class CompoundGroup:
    name: str


class Compound(SqlBaseClassComp):
    __tablename__ = "compound"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String, unique=True, nullable=False)

    smiles: Mapped[Optional[str]] = mapped_column(String)
    inchi: Mapped[Optional[str]] = mapped_column(String)
    inchkey: Mapped[Optional[str]] = mapped_column(String)
    formula: Mapped[Optional[str]] = mapped_column(String)

    rt: Mapped[Optional[float]] = mapped_column(Float)
    mz: Mapped[Optional[float]] = mapped_column(Float)
    ccs: Mapped[Optional[float]] = mapped_column(Float)

    # Relationship: Compound → IonPeak
    ions: Mapped[List["IonPeak"]] = relationship(
        back_populates="compound", cascade="all, delete-orphan"
    )


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


class FragmentPeak(SqlBaseClassComp, PeakBaseClass):
    __tablename__ = "fragment_peak"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    isotope_peak_id: Mapped[int] = mapped_column(ForeignKey("isotope_peak.id"), nullable=False)

    is_fragment: Mapped[bool] = mapped_column(Boolean, default=True)

    isotope_peak: Mapped["IsotopePeak"] = relationship(back_populates="fragments")


def test():
    # --- SQLITE in-memory engine ---
    engine = create_engine("sqlite:///:memory:", echo=False)

    # create tables
    SqlBaseClassComp.metadata.create_all(engine)

    # --- TEST INSERTION ---
    with Session(engine) as session:
        # Create a Compound
        cmp = Compound(
            name="TestCompound",
            formula="C10H20O",
            rt=5.5,
            mz=250.13,
            ccs=150.0
        )

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
            print(f"Compound: {c.name}, formula={c.formula}")
            for ion in c.ions:
                print(f"  IonPeak: mz={ion.mz}, adduct={ion.adduct}")
                for iso in ion.isotopes:
                    print(f"    IsotopePeak: mz={iso.mz}, order={iso.isotope_order}")
                    for frag in iso.fragments:
                        print(f"      FragmentPeak: mz={frag.mz}, intensity={frag.intensity}")


if __name__ == '__main__':
    test()
