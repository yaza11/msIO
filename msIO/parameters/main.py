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


rt = Parameter()
mz = Parameter()
mob = Parameter()
m = Parameter()
formula = Parameter()
adduct = Parameter()

