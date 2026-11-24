"""
Start the annotation GUI (msIO.gui.annotation): match measured features against a spectral library and annotate them.
"""
# noinspection PyUnresolvedReferences
import msIO.gui.executable  # noqa: F401  (must stay first)
import tkinter as tk

from msIO.gui.annotation import AnnotationGUI

if __name__ == "__main__":
    root = tk.Tk()
    app = AnnotationGUI(root)
    app.tk_root.mainloop()
