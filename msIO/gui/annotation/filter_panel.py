"""Filters of the results table, applied while typing."""

import tkinter as tk
from tkinter import ttk
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from LipidCalculator import CompoundDict

from msIO.gui.annotation.common import CollapsibleSection

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI


class FilterError(ValueError):
    """Invalid filter value, title and text can be shown to the user."""

    def __init__(self, title: str, text: str):
        super().__init__(text)
        self.title: str = title
        self.text: str = text


class FilterPanel(CollapsibleSection):
    """
    Range filters of the measured features (left) and filters of the matches (right). The filters are applied
    (AnnotationGUI.apply_filters) once typing pauses.
    """

    LIVE_FILTER_DELAY_MS: int = 300

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, "Filters (applied while typing)")
        self.app: "AnnotationGUI" = app
        self.live_filter_job: str = None
        self.last_values: tuple[str, ...] = None

        # each filter occupies a group of four grid columns: label, entry, "-", entry
        def add_label(row: int, group: int, text: str):
            ttk.Label(self.body, text=text).grid(row=row, column=4 * group, sticky="w", padx=(10, 4), pady=2)

        def add_range_filter(row: int, group: int, label: str) -> tuple[ttk.Entry, ttk.Entry]:
            add_label(row, group, label)
            tk_min = ttk.Entry(self.body, width=9)
            tk_min.grid(row=row, column=4 * group + 1, sticky="w")
            ttk.Label(self.body, text="-").grid(row=row, column=4 * group + 2, padx=2)
            tk_max = ttk.Entry(self.body, width=9)
            tk_max.grid(row=row, column=4 * group + 3, sticky="w")
            return tk_min, tk_max

        def add_filter(row: int, group: int, label: str) -> ttk.Entry:
            add_label(row, group, label)
            tk_entry = ttk.Entry(self.body)
            tk_entry.grid(row=row, column=4 * group + 1, columnspan=3, sticky="ew")
            return tk_entry

        # left: properties of the measured features, right: properties of the matches
        self.tk_min_f_id_meas, self.tk_max_f_id_meas = add_range_filter(0, 0, "Meas. feature ID")
        self.tk_min_rt, self.tk_max_rt = add_range_filter(1, 0, "RT (min)")
        self.tk_min_ccs, self.tk_max_ccs = add_range_filter(2, 0, "CCS (Å²)")
        self.tk_min_mz, self.tk_max_mz = add_range_filter(3, 0, "m/z")
        self.tk_min_intensity, self.tk_max_intensity = add_range_filter(4, 0, "Σ intensity")
        self.tk_min_ms2_score = add_filter(0, 1, "Min. MS2 score")
        self.tk_min_n_hits_ms2 = add_filter(1, 1, "Min. #MS2 peak matches")
        self.tk_subformula = add_filter(2, 1, "Formula contains")
        self.tk_name = add_filter(3, 1, "Name contains (a; b)")

        tk_buttons = ttk.Frame(self.body)
        tk_buttons.grid(row=4, column=4, columnspan=4, sticky="e")
        ttk.Button(tk_buttons, text="m/z from compound ...", command=app.on_click_mz_from_compound).pack(
            side="left", padx=(0, 5))
        ttk.Button(tk_buttons, text="Remove filters", command=app.on_click_remove_filters).pack(side="left")

        self.tk_entries: list[ttk.Entry] = [
            self.tk_min_f_id_meas, self.tk_max_f_id_meas, self.tk_name, self.tk_min_rt, self.tk_max_rt,
            self.tk_subformula, self.tk_min_ms2_score, self.tk_min_n_hits_ms2, self.tk_min_mz, self.tk_max_mz,
            self.tk_min_ccs, self.tk_max_ccs, self.tk_min_intensity, self.tk_max_intensity
        ]
        for tk_entry in self.tk_entries:
            tk_entry.bind("<KeyRelease>", self.on_edit)

    def get_values(self) -> tuple[str, ...]:
        return tuple(tk_entry.get().strip() for tk_entry in self.tk_entries)

    def on_edit(self, event=None):
        """Key release in a filter field: apply the filters once typing pauses."""
        # keys like Tab or arrows don't change the filters
        if self.get_values() == self.last_values:
            return
        # restart the timer on every key stroke, so filtering only happens once typing pauses
        if self.live_filter_job is not None:
            self.app.tk_root.after_cancel(self.live_filter_job)
        self.live_filter_job = self.app.tk_root.after(self.LIVE_FILTER_DELAY_MS, self._apply_live_filters)

    def _apply_live_filters(self):
        self.live_filter_job = None
        self.last_values = self.get_values()
        self.app.apply_filters(show_dialogs=False)

    def clear(self):
        for tk_entry in self.tk_entries:
            tk_entry.delete(0, tk.END)
        self.last_values = self.get_values()

    def set_mz_range(self, mz_min: float, mz_max: float):
        for tk_entry, value in ((self.tk_min_mz, mz_min), (self.tk_max_mz, mz_max)):
            tk_entry.delete(0, tk.END)
            tk_entry.insert(0, f"{value:.4f}")
        self.last_values = self.get_values()

    def compute_mask(self, table: pd.DataFrame) -> np.ndarray:
        """
        Rows of the matches table that pass all filters.

        Raises:
            FilterError: for an invalid filter value.
        """
        mask = np.ones(table.shape[0], dtype=bool)

        # numeric filters: (entry, column to compare, type, is minimum, error message)
        # features without a value (e.g. no CCS) are removed by filters on that column
        numeric_filters = [
            (self.tk_min_f_id_meas, table.feature_id_meas, int, True, 'Minimum ID must be a whole number.'),
            (self.tk_max_f_id_meas, table.feature_id_meas, int, False, 'Maximum ID must be a whole number.'),
            (self.tk_min_rt, table.RT, float, True, 'Minimum RT must be a number.'),
            (self.tk_max_rt, table.RT, float, False, 'Maximum RT must be a number.'),
            (self.tk_min_mz, table.mz_meas, float, True, 'Minimum m/z must be a number.'),
            (self.tk_max_mz, table.mz_meas, float, False, 'Maximum m/z must be a number.'),
            (self.tk_min_ms2_score, table.ms2_score, float, True, 'Minimum MS2 score must be a number.'),
            (self.tk_min_n_hits_ms2, table.n_hits_ms2, int, True, 'Minimum #MS2 peak matches must be a whole number.'),
            (self.tk_min_ccs, table.ccs_meas, float, True, 'Minimum CCS must be a number.'),
            (self.tk_max_ccs, table.ccs_meas, float, False, 'Maximum CCS must be a number.'),
            (self.tk_min_intensity, table.intensity, float, True, 'Minimum intensity must be a number.'),
            (self.tk_max_intensity, table.intensity, float, False, 'Maximum intensity must be a number.'),
        ]
        for tk_entry, values, value_type, is_minimum, error_message in numeric_filters:
            if not (text := tk_entry.get().strip()):
                continue
            try:
                limit = value_type(text)
            except ValueError:
                raise FilterError('Value Error', error_message)
            mask &= (values >= limit) if is_minimum else (values <= limit)

        # several name substrings can be separated by ";", all of them must be contained
        if name_substrings := self.tk_name.get().strip():
            for substring in name_substrings.split(';'):
                mask &= table.name.astype(str).str.contains(substring.strip(), case=False, na=False, regex=False)

        # formula must contain at least as many atoms of each element as the given subformula
        if subformula := self.tk_subformula.get().strip():
            try:
                composition_sub: dict[str, int] = CompoundDict(subformula).composition
                mask &= np.array([
                    all(CompoundDict(formula).composition.get(element, 0) >= count
                        for element, count in composition_sub.items())
                    for formula in table.formula
                ])
            except Exception as error:
                raise FilterError('Invalid Formula', f"Invalid formula: {error}")
        return np.asarray(mask, dtype=bool)
