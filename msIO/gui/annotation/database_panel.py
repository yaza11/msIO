"""Selection of the measured and the library database."""

import tkinter as tk
from tkinter import ttk
from typing import TYPE_CHECKING

from msIO.gui.annotation.common import CollapsibleSection

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI


class DatabasePanel(CollapsibleSection):
    """Paths of the measured and the library database, "Browse" loads them (see AnnotationGUI)."""

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, "Databases")
        self.app: "AnnotationGUI" = app
        self.tk_meas_path: tk.StringVar = tk.StringVar(app.tk_root)
        self.tk_lib_path: tk.StringVar = tk.StringVar(app.tk_root)

        ttk.Label(self.body, text="Measured database").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.body, textvariable=self.tk_meas_path).grid(row=0, column=1, padx=5, pady=2, sticky="ew")
        ttk.Button(self.body, text="Browse", command=app.on_browse_meas_db).grid(row=0, column=2)

        ttk.Label(self.body, text="Library database").grid(row=1, column=0, sticky="w")
        ttk.Entry(self.body, textvariable=self.tk_lib_path).grid(row=1, column=1, padx=5, pady=2, sticky="ew")
        ttk.Button(self.body, text="Browse", command=app.on_browse_lib_db).grid(row=1, column=2)
        self.body.columnconfigure(1, weight=1)
