import pint
from msIO.pint import ureg


class Parameter:
    default_name: str
    default_unit: pint.Unit
    names: list[str] = None
    default_symbol: str = None

    @classmethod
    def from_strings(cls, **ipt):
        ...


class PhysicalQuantity:
    name: str
    unit: str

    nominal_value: float
    uncertainty_absolute: float
    uncertainty_relative: float

    min_possible_value: float
    max_possible_value: float


# TODO: define common parameters
rt = Parameter()
mz = Parameter()
mob = Parameter()
m = Parameter()
formula = Parameter()
adduct = Parameter()

