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

    Dummy rows (one per measured feature, see matching.dummy_table) only hold feature properties: the feature filters
    apply to them, while a match filter is set, a dummy row is only kept if its feature has a match left. The maximum
    number of rows per feature is applied after filtering (AnnotationGUI.apply_filters), the dummy row comes first.
    """

    LIVE_FILTER_DELAY_MS: int = 300

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, "Filters (applied while typing)")
        self.app: "AnnotationGUI" = app
        self.live_filter_job: str = None
        self.last_values: tuple[str, ...] = None
        self.tk_add_dummies: tk.BooleanVar = tk.BooleanVar(app.tk_root, value=False)

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
        self.tk_min_i1_i0, self.tk_max_i1_i0 = add_range_filter(5, 0, "I1 / I0")
        self.tk_min_ms2_score = add_filter(0, 1, "Min. MS2 score")
        self.tk_min_n_hits_ms2 = add_filter(1, 1, "Min. #MS2 peak matches")
        self.tk_subformula = add_filter(2, 1, "Formula contains")
        self.tk_name = add_filter(3, 1, "Name contains (a; b)")
        self.tk_library = add_filter(4, 1, "Library contains (a | b)")

        # at most this many rows per feature (empty: all), optionally with a dummy row holding only the feature
        add_label(5, 1, "Max. rows per feature")
        self.tk_max_rows = ttk.Entry(self.body, width=9)
        self.tk_max_rows.grid(row=5, column=5, sticky="w")
        ttk.Checkbutton(self.body, text="Dummy row per feature", variable=self.tk_add_dummies,
                        command=lambda: app.apply_filters()).grid(row=5, column=6, columnspan=2, sticky="w",
                                                                  padx=(6, 0))

        tk_buttons = ttk.Frame(self.body)
        tk_buttons.grid(row=6, column=4, columnspan=4, sticky="e", pady=(2, 0))
        ttk.Button(tk_buttons, text="m/z from compound ...", command=app.on_click_mz_from_compound).pack(
            side="left", padx=(0, 5))
        ttk.Button(tk_buttons, text="Remove filters", command=app.on_click_remove_filters).pack(side="left")

        # filters of the match properties, see compute_mask
        self.tk_match_entries: list[ttk.Entry] = [
            self.tk_min_ms2_score, self.tk_min_n_hits_ms2, self.tk_subformula, self.tk_name, self.tk_library
        ]
        self.tk_entries: list[ttk.Entry] = [
            self.tk_min_f_id_meas, self.tk_max_f_id_meas, self.tk_min_rt, self.tk_max_rt, self.tk_min_ccs,
            self.tk_max_ccs, self.tk_min_mz, self.tk_max_mz, self.tk_min_intensity, self.tk_max_intensity,
            self.tk_min_i1_i0, self.tk_max_i1_i0, *self.tk_match_entries, self.tk_max_rows
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
        """Empty all filter fields (the dummy rows stay as they are)."""
        for tk_entry in self.tk_entries:
            tk_entry.delete(0, tk.END)
        self.last_values = self.get_values()

    def set_mz_range(self, mz_min: float, mz_max: float):
        for tk_entry, value in ((self.tk_min_mz, mz_min), (self.tk_max_mz, mz_max)):
            tk_entry.delete(0, tk.END)
            tk_entry.insert(0, f"{value:.4f}")
        self.last_values = self.get_values()

    def get_max_rows(self) -> int | None:
        """
        Maximum number of rows per measured feature, None for all.

        Raises:
            FilterError: if it is not a positive whole number.
        """
        if not (text := self.tk_max_rows.get().strip()):
            return None
        try:
            max_rows = int(text)
        except ValueError:
            max_rows = 0
        if max_rows < 1:
            raise FilterError('Value Error', 'Max. rows per feature must be a whole number of at least 1.')
        return max_rows

    @staticmethod
    def _apply_numeric_filters(table: pd.DataFrame, filters: list[tuple]) -> np.ndarray:
        """
        Rows within all limits, filters: (entry, column, type, is minimum, error message). Rows without a value (e.g.
        no CCS) are removed by filters on that column.
        """
        mask = np.ones(table.shape[0], dtype=bool)
        for tk_entry, column, value_type, is_minimum, error_message in filters:
            if not (text := tk_entry.get().strip()):
                continue
            try:
                limit = value_type(text)
            except ValueError:
                raise FilterError('Value Error', error_message)
            values = pd.to_numeric(table[column], errors='coerce')
            mask &= ((values >= limit) if is_minimum else (values <= limit)).to_numpy(dtype=bool, na_value=False)
        return mask

    def compute_mask(self, table: pd.DataFrame) -> np.ndarray:
        """
        Rows of the results table (matches and possibly dummy rows) that pass all filters.

        Raises:
            FilterError: for an invalid filter value.
        """
        feature_mask = self._apply_numeric_filters(table, [
            (self.tk_min_f_id_meas, 'feature_id_meas', int, True, 'Minimum ID must be a whole number.'),
            (self.tk_max_f_id_meas, 'feature_id_meas', int, False, 'Maximum ID must be a whole number.'),
            (self.tk_min_rt, 'RT', float, True, 'Minimum RT must be a number.'),
            (self.tk_max_rt, 'RT', float, False, 'Maximum RT must be a number.'),
            (self.tk_min_mz, 'mz_meas', float, True, 'Minimum m/z must be a number.'),
            (self.tk_max_mz, 'mz_meas', float, False, 'Maximum m/z must be a number.'),
            (self.tk_min_ccs, 'ccs_meas', float, True, 'Minimum CCS must be a number.'),
            (self.tk_max_ccs, 'ccs_meas', float, False, 'Maximum CCS must be a number.'),
            (self.tk_min_intensity, 'intensity', float, True, 'Minimum intensity must be a number.'),
            (self.tk_max_intensity, 'intensity', float, False, 'Maximum intensity must be a number.'),
            (self.tk_min_i1_i0, 'i1_i0', float, True, 'Minimum I1 / I0 must be a number.'),
            (self.tk_max_i1_i0, 'i1_i0', float, False, 'Maximum I1 / I0 must be a number.'),
        ])
        match_mask = self._apply_numeric_filters(table, [
            (self.tk_min_ms2_score, 'ms2_score', float, True, 'Minimum MS2 score must be a number.'),
            (self.tk_min_n_hits_ms2, 'n_hits_ms2', int, True, 'Minimum #MS2 peak matches must be a whole number.'),
        ])

        # several name substrings can be separated by ";", all of them must be contained
        if name_substrings := self.tk_name.get().strip():
            for substring in name_substrings.split(';'):
                match_mask &= table.name.astype(str).str.contains(
                    substring.strip(), case=False, na=False, regex=False).to_numpy(dtype=bool)

        # several libraries can be separated by "|", one of them must be contained in the file name of the library
        # or in the source of the library entry
        if library_substrings := [s.strip() for s in self.tk_library.get().split('|') if s.strip()]:
            texts: pd.Series = table.library.astype(str) + "\n" + table.source_library.astype(str)
            match_mask &= np.logical_or.reduce([
                texts.str.contains(substring, case=False, na=False, regex=False).to_numpy(dtype=bool)
                for substring in library_substrings
            ])

        # formula must contain at least as many atoms of each element as the given subformula
        if subformula := self.tk_subformula.get().strip():
            try:
                composition_sub: dict[str, int] = CompoundDict(subformula).composition
            except Exception as error:
                raise FilterError('Invalid Formula', f"Invalid formula: {error}")
            compositions: dict[str, dict[str, int]] = {}  # many matches share a formula

            def contains_subformula(formula: str) -> bool:
                if formula not in compositions:
                    try:
                        compositions[formula] = CompoundDict(formula).composition
                    except Exception:  # e.g. no formula
                        compositions[formula] = {}
                return all(compositions[formula].get(element, 0) >= count for element, count in composition_sub.items())

            match_mask &= np.array([contains_subformula(formula) for formula in table.formula.astype(str)], dtype=bool)

        is_dummy = table.is_dummy.to_numpy(dtype=bool) if 'is_dummy' in table.columns else np.zeros(len(table), bool)
        mask = feature_mask & match_mask & ~is_dummy
        if is_dummy.any():
            dummies = feature_mask & is_dummy
            if any(tk_entry.get().strip() for tk_entry in self.tk_match_entries):
                # dummy rows have no match properties: keep the ones of features with a match left
                dummies &= table.feature_id_meas.isin(set(table.feature_id_meas[mask])).to_numpy(dtype=bool)
            mask |= dummies
        return mask
