from importlib.util import find_spec as _find_spec

from ._version import __version__

# LipidCalculator is not on PyPI yet, so `pip install msIO` can't pull it in: fail with instructions instead of an
# ImportError deep inside msIO. find_spec only looks the package up without importing it (LipidCalculator imports msIO).
if _find_spec("LipidCalculator") is None:
    raise ImportError(
        "msIO needs LipidCalculator, install it with\n"
        "    pip install git+https://github.com/yaza11/LipidCalculator.git"
    )

from .list_of_ions.base import PeakList
from .list_of_ions.read_mgf import MgfImportManager
from .list_of_ions.read_msp import MSPReader
