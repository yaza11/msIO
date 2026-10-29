"""Filters of the results table, applied while typing."""

import re
import tkinter as tk
from tkinter import ttk
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from LipidCalculator import CompoundDict

from msIO.gui.annotation.common import FEATURE_RELATION_KINDS, CollapsibleSection

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI


class FilterError(ValueError):
    """Invalid filter value, title and text can be shown to the user."""

    def __init__(self, title: str, text: str):
        super().__init__(text)
        self.title: str = title
        self.text: str = text


_NUMBER: str = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
_RANGE = re.compile(rf'({_NUMBER})\s*-\s*({_NUMBER})?')
_TOLERANCE = re.compile(rf'({_NUMBER})\s*(?:±|\+-|\+/-)\s*({_NUMBER})\s*(ppm|%)?', re.IGNORECASE)
_BOUND = re.compile(rf'(<=|>=|<|>)\s*({_NUMBER})')
# (lower, upper, lower included, upper included)
Interval = tuple[float, float, bool, bool]
INTERVAL_HELP: str = ('Number fields take intervals, several separated by ";": 400-500; 512.3±0.01; 1302.4±5ppm; 7±1%; '
                      '>2; <=10; 3 (exactly)')


def parse_intervals(text: str) -> list[Interval]:
    """
    Intervals separated by ";":
        a-b          a to b (numbers may be negative: -5--2)
        a-           a and more
        c±t, c+-t    c - t to c + t, the tolerance t may be given in ppm or % of c (e.g. 1302.4±5ppm)
        >a, >=a, <b, <=b
        a            exactly a
    Empty text gives no intervals (no filter).

    Raises:
        ValueError: for a part that is none of these.
    """
    intervals: list[Interval] = []
    for part in text.split(';'):
        if not (part := part.strip()):
            continue
        if match := _TOLERANCE.fullmatch(part):
            center, tolerance = float(match.group(1)), abs(float(match.group(2)))
            if unit := (match.group(3) or '').lower():
                tolerance *= abs(center) * (1e-6 if unit == 'ppm' else 1e-2)
            intervals.append((center - tolerance, center + tolerance, True, True))
        elif match := _RANGE.fullmatch(part):
            lower = float(match.group(1))
            upper = np.inf if match.group(2) is None else float(match.group(2))
            if upper < lower:
                raise ValueError(f'"{part}": the upper limit is below the lower one')
            intervals.append((lower, upper, True, True))
        elif match := _BOUND.fullmatch(part):
            operator, value = match.group(1), float(match.group(2))
            if operator.startswith('>'):
                intervals.append((value, np.inf, operator == '>=', False))
            else:
                intervals.append((-np.inf, value, False, operator == '<='))
        else:
            try:
                value = float(part)
            except ValueError:
                raise ValueError(f'"{part}" is no interval (e.g. 400-500, 512.3±0.01, >2 or 3)') from None
            intervals.append((value, value, True, True))
    return intervals


def in_intervals(values: pd.Series, intervals: list[Interval]) -> np.ndarray:
    """Values in any of the intervals, missing values are in none."""
    numbers: np.ndarray = pd.to_numeric(values, errors='coerce').to_numpy(dtype=float, na_value=np.nan)
    mask = np.zeros(len(numbers), dtype=bool)
    with np.errstate(invalid='ignore'):  # nan
        for lower, upper, lower_included, upper_included in intervals:
            above = numbers >= lower if lower_included else numbers > lower
            below = numbers <= upper if upper_included else numbers < upper
            mask |= above & below
    return mask


def split_terms(text: str, separators: str) -> list[str]:
    return [term.strip() for term in re.split(f'[{re.escape(separators)}]', text) if term.strip()]


class FilterPanel(CollapsibleSection):
    """
    Filters of the measured features (left), of the match properties (middle) and of the match texts (right). The
    filters are applied (AnnotationGUI.apply_filters) once typing pauses. Number fields take intervals (see
    parse_intervals), the "Min." fields a lower limit.

    Dummy rows (one per measured feature, see matching.dummy_table) only hold feature properties: the feature filters
    apply to them, while a match filter is set (blacklists don't count), a dummy row is only kept if its feature has a
    match left. The maximum number of rows per feature is applied after filtering (AnnotationGUI.apply_filters), the
    dummy row comes first.
    """

    LIVE_FILTER_DELAY_MS: int = 300

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, "Filters (applied while typing)")
        self.app: "AnnotationGUI" = app
        self.live_filter_job: str = None
        self.last_values: tuple[str, ...] = None
        self.tk_add_dummies: tk.BooleanVar = tk.BooleanVar(app.tk_root, value=False)
        self.tk_hide_related: tk.BooleanVar = tk.BooleanVar(app.tk_root, value=False)

        # each filter occupies a group of two grid columns: label, entry
        def add_filter(row: int, group: int, label: str, width: int = 16) -> ttk.Entry:
            ttk.Label(self.body, text=label).grid(row=row, column=2 * group, sticky="w", padx=(10, 4), pady=2)
            tk_entry = ttk.Entry(self.body, width=width)
            tk_entry.grid(row=row, column=2 * group + 1, sticky="w")
            return tk_entry

        # properties of the measured features: intervals
        self.tk_f_id_meas = add_filter(0, 0, "Meas. feature ID")
        self.tk_rt = add_filter(1, 0, "RT (min)")
        self.tk_ccs = add_filter(2, 0, "CCS (Å²)")
        self.tk_mz = add_filter(3, 0, "m/z")
        self.tk_intensity = add_filter(4, 0, "Σ intensity")
        self.tk_i1_i0 = add_filter(5, 0, "I1 / I0")
        # numbers of the matches
        self.tk_min_ms2_score = add_filter(0, 1, "Min. MS2 score", width=10)
        self.tk_min_n_hits_ms2 = add_filter(1, 1, "Min. #MS2 peak matches", width=10)
        self.tk_dmz_mda = add_filter(2, 1, "|Δm/z| (mDa)", width=10)
        self.tk_dmz_ppm = add_filter(3, 1, "|Δm/z| (ppm)", width=10)
        self.tk_dccs = add_filter(4, 1, "ΔCCS (%)", width=10)
        self.tk_xlogp = add_filter(5, 1, "logP", width=10)
        # texts of the matches
        self.tk_subformula = add_filter(0, 2, "Formula contains", width=18)
        self.tk_name = add_filter(1, 2, "Name contains (a; b)", width=18)
        self.tk_name_excluded = add_filter(2, 2, "Name excludes (a; b)", width=18)
        self.tk_library = add_filter(3, 2, "Library contains (a | b)", width=18)
        self.tk_library_excluded = add_filter(4, 2, "Library excludes (a; b)", width=18)

        # at most this many rows per feature (empty: all), optionally with a dummy row holding only the feature
        tk_rows = ttk.Frame(self.body)
        tk_rows.grid(row=5, column=4, columnspan=2, sticky="w", padx=(10, 0))
        ttk.Label(tk_rows, text="Max. rows per feature").pack(side="left", padx=(0, 4))
        self.tk_max_rows = ttk.Entry(tk_rows, width=6)
        self.tk_max_rows.pack(side="left")
        ttk.Checkbutton(tk_rows, text="Dummy row per feature", variable=self.tk_add_dummies,
                        command=lambda: app.apply_filters()).pack(side="left", padx=(8, 0))

        ttk.Label(self.body, text=INTERVAL_HELP, foreground="gray").grid(row=6, column=0, columnspan=4, sticky="w",
                                                                         padx=(10, 0), pady=(2, 0))
        tk_buttons = ttk.Frame(self.body)
        tk_buttons.grid(row=7, column=0, columnspan=6, sticky="ew", pady=(2, 0))
        ttk.Checkbutton(tk_buttons, text="Hide features marked as fragment, duplicate, isotope peak or adduct",
                        variable=self.tk_hide_related, command=lambda: app.apply_filters()).pack(
            side="left", padx=(10, 0))
        ttk.Button(tk_buttons, text="Remove filters", command=app.on_click_remove_filters).pack(side="right")
        ttk.Button(tk_buttons, text="m/z from compound ...", command=app.on_click_mz_from_compound).pack(
            side="right", padx=5)

        # (entry, column of the results table, name in error messages)
        self.feature_intervals: list[tuple[ttk.Entry, str, str]] = [
            (self.tk_f_id_meas, 'feature_id_meas', 'Meas. feature ID'),
            (self.tk_rt, 'RT', 'RT'),
            (self.tk_ccs, 'ccs_meas', 'CCS'),
            (self.tk_mz, 'mz_meas', 'm/z'),
            (self.tk_intensity, 'intensity', 'Σ intensity'),
            (self.tk_i1_i0, 'i1_i0', 'I1 / I0'),
        ]
        self.match_intervals: list[tuple[ttk.Entry, str, str]] = [
            (self.tk_dmz_mda, 'dmz_mda', '|Δm/z| (mDa)'),
            (self.tk_dmz_ppm, 'dmz_ppm', '|Δm/z| (ppm)'),
            (self.tk_dccs, 'dccs_percent', 'ΔCCS'),
            (self.tk_xlogp, 'xlogp', 'logP'),
        ]
        # filters that require a match (see compute_mask)
        self.tk_match_entries: list[ttk.Entry] = [
            self.tk_min_ms2_score, self.tk_min_n_hits_ms2, *(entry for entry, _, _ in self.match_intervals),
            self.tk_subformula, self.tk_name, self.tk_library
        ]
        self.tk_entries: list[ttk.Entry] = [
            *(entry for entry, _, _ in self.feature_intervals), *self.tk_match_entries, self.tk_name_excluded,
            self.tk_library_excluded, self.tk_max_rows
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
        """Empty all filter fields (the check boxes stay as they are)."""
        for tk_entry in self.tk_entries:
            tk_entry.delete(0, tk.END)
        self.last_values = self.get_values()

    def set_mz_text(self, text: str):
        """Replace the m/z filter (intervals, see parse_intervals)."""
        self.tk_mz.delete(0, tk.END)
        self.tk_mz.insert(0, text)
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
    def _apply_interval_filters(table: pd.DataFrame, filters: list[tuple[ttk.Entry, str, str]]) -> np.ndarray:
        """Rows within the intervals of all filters, rows without a value (e.g. no CCS) are removed by its filter."""
        mask = np.ones(table.shape[0], dtype=bool)
        for tk_entry, column, name in filters:
            try:
                intervals: list[Interval] = parse_intervals(tk_entry.get())
            except ValueError as error:
                raise FilterError('Value Error', f'{name}: {error}')
            if intervals:
                mask &= in_intervals(table[column], intervals)
        return mask

    @staticmethod
    def _apply_minimum_filters(table: pd.DataFrame, filters: list[tuple]) -> np.ndarray:
        """Rows with at least the minimum, filters: (entry, column, type, error message)."""
        mask = np.ones(table.shape[0], dtype=bool)
        for tk_entry, column, value_type, error_message in filters:
            if not (text := tk_entry.get().strip()):
                continue
            try:
                limit = value_type(text)
            except ValueError:
                raise FilterError('Value Error', error_message)
            values = pd.to_numeric(table[column], errors='coerce')
            mask &= (values >= limit).to_numpy(dtype=bool, na_value=False)
        return mask

    @staticmethod
    def _contains_any(texts: pd.Series, substrings: list[str]) -> np.ndarray:
        return np.logical_or.reduce([
            texts.str.contains(substring, case=False, na=False, regex=False).to_numpy(dtype=bool)
            for substring in substrings
        ])

    def compute_mask(self, table: pd.DataFrame) -> np.ndarray:
        """
        Rows of the results table (matches and possibly dummy rows) that pass all filters.

        Raises:
            FilterError: for an invalid filter value.
        """
        feature_mask = self._apply_interval_filters(table, self.feature_intervals)
        if self.tk_hide_related.get():
            related: set[int] = self.app.annotations.annotated_feature_ids(kinds=FEATURE_RELATION_KINDS)
            feature_mask &= ~table.feature_id_meas.isin(related).to_numpy(dtype=bool)

        match_mask = self._apply_interval_filters(table, self.match_intervals)
        match_mask &= self._apply_minimum_filters(table, [
            (self.tk_min_ms2_score, 'ms2_score', float, 'Minimum MS2 score must be a number.'),
            (self.tk_min_n_hits_ms2, 'n_hits_ms2', int, 'Minimum #MS2 peak matches must be a whole number.'),
        ])

        # several name substrings can be separated by ";", all of them must be contained, none of the excluded ones
        names: pd.Series = table.name.astype(str)
        for substring in split_terms(self.tk_name.get(), ';'):
            match_mask &= names.str.contains(substring, case=False, na=False, regex=False).to_numpy(dtype=bool)
        if excluded := split_terms(self.tk_name_excluded.get(), ';|'):
            match_mask &= ~self._contains_any(names, excluded)

        # several libraries can be separated by "|", one of them must be contained in the file name of the library
        # or in the source of the library entry, none of the excluded ones
        libraries: pd.Series = table.library.astype(str) + "\n" + table.source_library.astype(str)
        if library_substrings := split_terms(self.tk_library.get(), '|'):
            match_mask &= self._contains_any(libraries, library_substrings)
        if excluded := split_terms(self.tk_library_excluded.get(), ';|'):
            match_mask &= ~self._contains_any(libraries, excluded)

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
