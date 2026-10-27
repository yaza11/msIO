"""Parameters of the library search."""

import tkinter as tk
from tkinter import ttk
from typing import TYPE_CHECKING

from msIO.gui.annotation.common import CollapsibleSection

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI


class MatchingPanel(CollapsibleSection):
    """
    Tolerances, MS2 score metric and the "Find Matches" button. The precursor tolerance is the stricter of the
    ppm and mDa tolerance, the tolerances are also used by the overview and the scatter plot.
    """

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, "Matching parameters")
        self.app: "AnnotationGUI" = app
        self.tk_max_dmz_ppm: tk.DoubleVar = tk.DoubleVar(app.tk_root, value=10)
        self.tk_max_dmz_mda: tk.DoubleVar = tk.DoubleVar(app.tk_root, value=3)
        self.tk_max_dmz_ms2_mda: tk.DoubleVar = tk.DoubleVar(app.tk_root, value=10)
        self.tk_min_ms2_score: tk.DoubleVar = tk.DoubleVar(app.tk_root, value=0)
        self.tk_require_ms2: tk.BooleanVar = tk.BooleanVar(app.tk_root, value=False)

        def add_entry(row: int, column: int, label: str, variable: tk.Variable):
            ttk.Label(self.body, text=label).grid(row=row, column=column, sticky="w", padx=(10, 4), pady=2)
            ttk.Entry(self.body, textvariable=variable, width=8).grid(row=row, column=column + 1, sticky="w")

        add_entry(0, 0, "Precursor m/z tol. (ppm)", self.tk_max_dmz_ppm)
        add_entry(0, 2, "and (mDa)", self.tk_max_dmz_mda)
        add_entry(1, 0, "MS2 fragment m/z tol. (mDa)", self.tk_max_dmz_ms2_mda)
        add_entry(1, 2, "Min. MS2 score", self.tk_min_ms2_score)

        ttk.Label(self.body, text="MS2 score metric").grid(row=0, column=4, sticky="w", padx=(10, 4), pady=2)
        self.tk_metric: ttk.Combobox = ttk.Combobox(
            self.body, state="readonly", values=['cosine_fwd', 'cosine_bwd', 'cosine_sim', 'modified_cosine_greedy'],
            width=22
        )
        self.tk_metric.set('modified_cosine_greedy')
        self.tk_metric.grid(row=0, column=5, sticky="w")
        ttk.Checkbutton(self.body, text="Require MS2 match", variable=self.tk_require_ms2).grid(
            row=1, column=4, columnspan=2, sticky="w", padx=(10, 0)
        )

        self.tk_find_button = ttk.Button(self.body, text="Find Matches", command=app.on_click_find_matches)
        self.tk_find_button.grid(row=0, column=6, rowspan=2, sticky="nse", padx=(20, 0))
        self.body.columnconfigure(6, weight=1)

    def get_search_parameters(self) -> dict:
        """
        Keyword arguments for Library.find_matches (without mzs and spectra).

        Raises:
            ValueError: if a parameter is invalid (the message can be shown to the user).
        """
        try:
            parameters = dict(
                max_dmz_ppm=self.tk_max_dmz_ppm.get(),
                max_dmz_da=self.tk_max_dmz_mda.get() * 1e-3,
                max_ms2_dmz_da=self.tk_max_dmz_ms2_mda.get() * 1e-3,
                min_ms2_score=self.tk_min_ms2_score.get(),
                metric=self.tk_metric.get(),
                require_ms2=self.tk_require_ms2.get(),
            )
        except tk.TclError:
            raise ValueError("Matching parameters must be numbers.")
        if parameters['max_dmz_ppm'] <= 0 or parameters['max_dmz_da'] <= 0 or parameters['max_ms2_dmz_da'] <= 0:
            raise ValueError("Tolerances must be greater than 0.")
        if not 0 <= parameters['min_ms2_score'] <= 1:
            raise ValueError("Please enter a MS2 score between 0 and 1.")
        return parameters

    def precursor_tolerances(self) -> tuple[float, float]:
        """(ppm, Da), raises tk.TclError for invalid input."""
        return self.tk_max_dmz_ppm.get(), self.tk_max_dmz_mda.get() * 1e-3

    def precursor_tolerance_da(self, mz: float) -> float:
        """Precursor m/z tolerance in Da at mz: the stricter of the ppm and mDa tolerance (as in the search)."""
        max_dmz_ppm, max_dmz_da = self.precursor_tolerances()
        return min(max_dmz_da, max_dmz_ppm * 1e-6 * mz)

    def ms2_tolerance_da(self) -> float:
        """MS2 fragment tolerance in Da, raises tk.TclError for invalid input."""
        return self.tk_max_dmz_ms2_mda.get() * 1e-3

    def set_searching(self, searching: bool):
        """Disable "Find Matches" while the library is searched (or indexed)."""
        self.tk_find_button.config(state="disabled" if searching else "normal")

    def is_searching(self) -> bool:
        return str(self.tk_find_button.cget("state")) == "disabled"
