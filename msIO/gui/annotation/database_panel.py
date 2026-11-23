"""Selection of the measured database and the library databases."""

import tkinter as tk
from tkinter import ttk
from typing import TYPE_CHECKING

from msIO.gui.annotation.common import CollapsibleSection, add_tooltip

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI


class DatabasePanel(CollapsibleSection):
    """
    Path of the measured database and paths of the libraries (separated by ";"), see AnnotationGUI:
    "Browse" selects the measured database / replaces the libraries, "Add" adds libraries, Return in the library field
    loads the typed paths.
    """

    SEPARATOR: str = ";"

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, "Databases", tooltip="The measured database (Ctrl+O) and the library databases "
                                                      "(Ctrl+Shift+O: add) to search.")
        self.app: "AnnotationGUI" = app
        self.tk_meas_path: tk.StringVar = tk.StringVar(app.tk_root)
        self.tk_lib_paths: tk.StringVar = tk.StringVar(app.tk_root)

        tk_meas_label = ttk.Label(self.body, text="Measured database")
        tk_meas_label.grid(row=0, column=0, sticky="w")
        add_tooltip("SQLite file of the measurement (e.g. from the MZmine/MetaboScape converter). Ctrl+O: open.",
                    tk_meas_label)
        ttk.Entry(self.body, textvariable=self.tk_meas_path).grid(row=0, column=1, padx=5, pady=2, sticky="ew")
        ttk.Button(self.body, text="Browse", command=app.on_browse_meas_db).grid(row=0, column=2)

        tk_lib_label = ttk.Label(self.body, text="Libraries")
        tk_lib_label.grid(row=1, column=0, sticky="w")
        add_tooltip("Library databases, several separated by \";\" (Enter loads typed paths). Browse: replace, Add / "
                    "Ctrl+Shift+O: add libraries.", tk_lib_label)
        tk_lib_entry = ttk.Entry(self.body, textvariable=self.tk_lib_paths)
        tk_lib_entry.grid(row=1, column=1, padx=5, pady=2, sticky="ew")
        tk_lib_entry.bind("<Return>", lambda event: app.load_libraries(self.get_lib_paths()))
        ttk.Button(self.body, text="Browse", command=app.on_browse_lib_db).grid(row=1, column=2)
        ttk.Button(self.body, text="Add", command=app.on_add_lib_db).grid(row=1, column=3, padx=(3, 0))
        self.body.columnconfigure(1, weight=1)

    def get_lib_paths(self) -> list[str]:
        return [path.strip() for path in self.tk_lib_paths.get().split(self.SEPARATOR) if path.strip()]

    def set_lib_paths(self, paths: list[str]):
        self.tk_lib_paths.set(f"{self.SEPARATOR} ".join(paths))
