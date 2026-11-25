"""
Give windowed executables a stdout / stderr.

Windowed PyInstaller builds (--windowed) and pythonw have no console: sys.stdout / sys.stderr are None, so print(),
tqdm, logging and tracebacks of exceptions in tk callbacks would fail ("'NoneType' object has no attribute 'write'")
or vanish. Both streams are redirected to a log file in the user's home directory (devnull if that is not writable).

Used in two ways:
    - as PyInstaller runtime hook (--runtime-hook msIO/gui/_gui_log.py), runs before any other import in the exe
    - imported first by every GUI entry module, covers pythonw and running the modules without PyInstaller

Only uses the standard library and is idempotent, so running it twice does no harm.
"""
import os
import sys

LOG_FILE = os.path.join(os.path.expanduser("~"), "msIO_gui.log")


def _open_log():
    try:
        return open(LOG_FILE, "a", buffering=1, encoding="utf-8", errors="replace")
    except OSError:
        return open(os.devnull, "w")


if sys.stdout is None or sys.stderr is None:
    _log = _open_log()
    if sys.stdout is None:
        sys.stdout = _log
    if sys.stderr is None:
        sys.stderr = _log
