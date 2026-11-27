"""Temporary launcher for the Binder desktop.

master imports msIO.gui.executable, but only msIO/gui/_gui_log.py exists (commit 149a250). This maps the old name to
the new module, so the GUI starts either way. Remove this file once the imports point to msIO.gui._gui_log.
"""
import runpy
import sys

try:
    import msIO.gui.executable  # noqa: F401
except ModuleNotFoundError:
    import msIO.gui._gui_log as _gui_log
    sys.modules["msIO.gui.executable"] = _gui_log

runpy.run_module("msIO.gui.main", run_name="__main__")