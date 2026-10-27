"""
Start the annotation GUI (msIO.gui.annotation): match measured features against a spectral library and annotate them.
"""

import tkinter as tk

from msIO.gui.annotation import AnnotationGUI

if __name__ == "__main__":
    root = tk.Tk()
    app = AnnotationGUI(root)
    app.tk_root.mainloop()
