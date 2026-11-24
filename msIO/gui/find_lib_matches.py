"""
Start the annotation GUI (msIO.gui.annotation): match measured features against a spectral library and annotate them.
"""
import os
import sys
if sys.stdout is None:
    sys.stdout = open(os.path.join(os.path.expanduser("~"), "msio_converter.log"), "a", buffering=1)
if sys.stderr is None:
    sys.stderr = open(os.path.join(os.path.expanduser("~"), "msio_converter_errors.log"), "a", buffering=1)

import msIO.gui.executable

import tkinter as tk

from msIO.gui.annotation import AnnotationGUI

if __name__ == "__main__":
    root = tk.Tk()
    app = AnnotationGUI(root)
    app.tk_root.mainloop()
