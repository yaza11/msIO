from sqlalchemy import Float, String, Boolean, Integer
from sqlalchemy.orm import Mapped, mapped_column

from msIO.features.base import FeatureBaseClass


class PeakBaseClass(FeatureBaseClass):
    id: Mapped[int] = mapped_column(primary_key=True)
    mz: Mapped[float] = mapped_column(Float)
    rt: Mapped[float | None] = mapped_column(Float, nullable=True)
    ccs: Mapped[float | None] = mapped_column(Float, nullable=True)

    intensity: Mapped[float | None] = mapped_column(Float, nullable=True)
    height: Mapped[float | None] = mapped_column(Float, nullable=True)
    area: Mapped[float | None] = mapped_column(Float, nullable=True)

    fwhm: Mapped[float | None] = mapped_column(Float, nullable=True)
    snr: Mapped[float | None] = mapped_column(Float, nullable=True)

    M: Mapped[float | None] = mapped_column(Float, nullable=True)
    adduct: Mapped[str | None] = mapped_column(String, nullable=True)

    is_ion: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    is_isotope: Mapped[bool | None] = mapped_column(Boolean, nullable=True)

    is_fragment: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
