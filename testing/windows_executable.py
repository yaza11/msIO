import os
import sys
if sys.stdout is None or sys.stderr is None:
    _log = open(os.path.join(os.path.expanduser("~"), "msio_converter.log"), "a", buffering=1)
    if sys.stdout is None:
        sys.stdout = _log
    if sys.stderr is None:
        sys.stderr = _log