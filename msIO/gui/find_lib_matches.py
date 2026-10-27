"""
GUI for annotating measured features with matches from a spectral library.

Workflow:
    1. Select a measured database (FeatureManagerDB) and a library database (Library).
       The scatter plot (m/z vs RT) of the measurement is drawn right away and the MS2
       spectra of the measurement are loaded in the background.
    2. Set the matching parameters and click "Find Matches".
    3. Narrow down the results table with the filters (applied while typing). Clicking a
       feature in the scatter plot selects its best match in the table. In the scatter plot, features with specific
       MS2 fragments or neutral losses can be highlighted and pairs of features with
       specific mass differences (e.g. adduct exchanges) can be connected.
    4. The overview panel shows the selected match (table row), a single feature (click in
       the scatter plot or feature ID field) or two compared features (shift-click two
       features or enter two IDs): MS2 spectra (peaks can be annotated by clicking them),
       MS1 spectrum, intensities in the samples and the structure (see feature_overview).
    5. Collect annotations in a local list: library matches ("Add library match"),
       compounds entered by hand or built with the LipidCalculator IPL builder ("Add manual
       annotation") and in-source fragments of other features ("Mark in-source fragment").
       "Manage annotations" opens a window to review new and existing annotations and to
       write the changes to the measured database.

Slow work (loading MS2 spectra, building the match overview figure) runs in worker
threads. Worker threads never touch tk widgets; they hand their results back to the
main thread with ``tk_root.after``.
"""

import numpy as np
import os
import pandas as pd
import re
import threading
import tkinter as tk
from tkinter import messagebox, filedialog, simpledialog, ttk
from typing import Callable

from LipidCalculator import Adduct, CompoundDict
from LipidCalculator.compound_creation.formula_parser import parse_equation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from msIO.feature_managers.db import Library, FeatureManagerDB, lib_match_result_to_table
from msIO.gui.feature_overview import (
    OverviewData, OverviewSettings, PlottedSpectrum, SpectrumKey, build_overview_figure, gather_feature,
    gather_library, peak_key
)

try:
    from LipidCalculator.gui.common import ADDUCTS
except ImportError:  # LipidCalculator without GUI
    ADDUCTS: list[str] = ['[M+H]+', '[M+NH4]+', '[M+Na]+', '[M+K]+', '[M-H]-', '[M+Cl]-', '[M+HCOO]-']

SQLITE_FILE_TYPES: list[tuple[str, str]] = [("SQLite Files", "*.sqlite *.db *.sql"), ("All files", "*.*")]

# (kind, measured feature id, reference) identifies a manual annotation, kind is one of the
# ManualAnnotationKind of FeatureManagerDB, the reference is the library feature id ("library"),
# the parent feature id ("in_source_fragment") or "name|adduct" ("manual")
AnnotationKey = tuple[str, int, str]
ANNOTATION_KIND_LABELS: dict[str, str] = {
    "library": "library match", "manual": "manual", "in_source_fragment": "in-source fragment"
}

# colors of the two features compared by shift-click, used in the scatter plot and the MS2 comparison
COMPARISON_COLORS: tuple[str, str] = ("C2", "C1")
# matplotlib key names of + and - (main keyboard and numpad) for zooming the scatter plot
ZOOM_IN_KEYS: set[str] = {"+", "plus", "add", "kp_add"}
ZOOM_OUT_KEYS: set[str] = {"-", "minus", "subtract", "kp_subtract"}
# colors of the lines connecting features with marked mass differences (one per entry)
MASS_DIFFERENCE_COLORS: list[str] = ["C3", "C4", "C5", "C6", "C8", "C9"]


def parse_masses(text: str) -> list[tuple[str, float]]:
    """
    Parse ";"-separated m/z values or formulas into (entry, mass) pairs.

    Formulas are converted with LipidCalculator, a trailing "+" or "-" gives the mass of the ion
    (e.g. "C5H15NO4P+" for the phosphocholine fragment).

    Raises:
        ValueError: if an entry is neither a number nor a valid formula, or its mass is not positive.
    """
    masses: list[tuple[str, float]] = []
    for entry in text.split(';'):
        if not (entry := entry.strip()):
            continue
        try:
            mass = float(entry)
        except ValueError:
            try:
                mass = CompoundDict(entry).mass
            except Exception:
                raise ValueError(f'"{entry}" is neither a number nor a valid formula.')
        if mass <= 0:
            raise ValueError(f'"{entry}" does not have a positive mass.')
        masses.append((entry, mass))
    return masses


def parse_mass_differences(text: str) -> list[tuple[str, Callable[[np.ndarray], np.ndarray]]]:
    """
    Parse ";"-separated mass differences into (entry, function) pairs. The function maps m/z values
    of features to the m/z values of their expected partners.

    An entry is a mass difference (number, formula or equation like "C[13] - C[12]" or "H2O + NH3", see
    LipidCalculator's parse_equation) or an adduct exchange like "[M+H]+>[M+NH4]+" (the partner of an
    [M+H]+ ion is the [M+NH4]+ ion of the same molecule, this also works for different charge states).

    Raises:
        ValueError: if an entry can't be parsed.
    """
    differences: list[tuple[str, Callable[[np.ndarray], np.ndarray]]] = []
    for entry in text.split(';'):
        if not (entry := entry.strip()):
            continue
        if '>' in entry:
            try:
                adduct_from, adduct_to = (Adduct(adduct.strip()) for adduct in entry.split('>'))
            except Exception:
                raise ValueError(f'"{entry}" is not a valid adduct exchange, use e.g. [M+H]+>[M+NH4]+.')
            differences.append((entry, lambda mzs, a=adduct_from, b=adduct_to: np.array(
                [b.mass_to_mz(a.mz_to_mass(mz)) for mz in mzs]
            )))
        else:
            try:
                mass = float(entry)
            except ValueError:
                try:
                    mass = parse_equation(entry).mass
                except Exception:
                    raise ValueError(f'"{entry}" is neither a number, a formula nor an equation like C[13] - C[12].')
            if mass == 0:
                raise ValueError(f'"{entry}" has no mass difference.')
            differences.append((entry, lambda mzs, d=mass: mzs + d))
    return differences


def open_ipl_builder(parent: tk.Misc, send_label: str, on_send: Callable[[object, dict, tk.Toplevel], None]):
    """
    Open the IPL builder of LipidCalculator in a window. Its button send_label calls on_send(frame, payload, window),
    the payload holds the neutral formula and the adduct selected in the table of the builder (if any).
    """
    try:
        from LipidCalculator.gui.ipl_gui import SingleIPLFrame
    except Exception as error:
        messagebox.showerror("IPL builder", f"The IPL builder of LipidCalculator is not available:\n{error}",
                             parent=parent)
        return
    tk_window = tk.Toplevel(parent)
    tk_window.title(f"IPL builder (select an adduct in the table, then \"{send_label}\")")
    tk_window.geometry("1100x750")
    frames: list = []  # the frame doesn't exist yet when the callback is defined
    frame = SingleIPLFrame(tk_window, send_targets={send_label: lambda payload: on_send(frames[0], payload, tk_window)})
    frames.append(frame)
    frame.pack(fill="both", expand=True)


def find_feature_pairs(
        f_ids: np.ndarray,
        rts: np.ndarray,
        mzs: np.ndarray,
        expected_partner_mz: Callable[[np.ndarray], np.ndarray],
        max_dmz_ppm: float,
        max_dmz_da: float,
        rt_window: float
) -> list[tuple[int, int]]:
    """
    Pairs (feature, partner) whose partner has the expected m/z (stricter of both tolerances) and a
    retention time within rt_window of the feature.
    """
    expected: np.ndarray = expected_partner_mz(mzs)
    tolerances: np.ndarray = np.minimum(max_dmz_da, max_dmz_ppm * 1e-6 * expected)
    order: np.ndarray = np.argsort(mzs)
    mzs_sorted: np.ndarray = mzs[order]
    idcs_left = np.searchsorted(mzs_sorted, expected - tolerances, side='left')
    idcs_right = np.searchsorted(mzs_sorted, expected + tolerances, side='right')

    pairs: list[tuple[int, int]] = []
    for i, (left, right) in enumerate(zip(idcs_left, idcs_right)):
        for j in order[left:right]:
            if j != i and abs(rts[j] - rts[i]) <= rt_window:
                pairs.append((int(f_ids[i]), int(f_ids[j])))
    return pairs


class AnnotationGUI:
    """
    Main window for matching measured features against a library and annotating them.

    Naming conventions:
        tk_*        tk widgets and tk variables
        on_*        callbacks of widgets and events
        _*worker    functions running in a worker thread
    """

    def __init__(self, root: tk.Tk):
        self.tk_root: tk.Tk = root
        self.tk_root.title("Annotation GUI")
        # 80 % of the screen, the plots need the space
        width: int = int(self.tk_root.winfo_screenwidth() * 0.8)
        height: int = int(self.tk_root.winfo_screenheight() * 0.8)
        self.tk_root.geometry(f"{width}x{height}+10+10")
        self.tk_root.protocol("WM_DELETE_WINDOW", self.on_close)
        # tk variables are created with root as master: without it they belong to the first Tk instance of the process,
        # which is not this window if e.g. an IPython console created one before, and the entries would stay empty

        # database selection
        self.tk_meas_path: tk.StringVar = tk.StringVar(self.tk_root)
        self.tk_lib_path: tk.StringVar = tk.StringVar(self.tk_root)
        self.meas_db: FeatureManagerDB = None
        self.lib_db: Library = None
        # incremented for every loaded measured database, so results of workers for an older database are dropped
        self.meas_db_load_number: int = 0

        # MS2 spectra of the measurement, loaded in a worker thread
        self.ms2_spectra: dict = None
        self.ms2_loading: bool = False
        self.ms2_loading_error: str = None

        # matching parameters
        self.tk_max_dmz_ppm: tk.DoubleVar = tk.DoubleVar(self.tk_root, value=10)
        self.tk_max_dmz_mda: tk.DoubleVar = tk.DoubleVar(self.tk_root, value=3)
        self.tk_max_dmz_ms2_mda: tk.DoubleVar = tk.DoubleVar(self.tk_root, value=10)
        self.tk_min_ms2_score: tk.DoubleVar = tk.DoubleVar(self.tk_root, value=0)
        self.tk_require_ms2: tk.BooleanVar = tk.BooleanVar(self.tk_root, value=False)

        # matching results: raw output of Library.find_matches and as table (all and filtered)
        self.matches: dict[int, list[dict]] = {}
        self.matches_table: pd.DataFrame = pd.DataFrame()
        self.filtered_matches_table: pd.DataFrame = pd.DataFrame()
        # sorting of the results table: (column, descending) per clicked header, the first one is the primary key,
        # kept when the filters change; empty: default order
        self.sort_keys: list[tuple[str, bool]] = []

        # filters are applied while typing, after a short pause
        self.live_filter_delay_ms: int = 300
        self.live_filter_job: str = None
        self.last_filter_values: tuple[str, ...] = None

        # manual annotations: changes are collected locally and only written to the measured database on request
        self.saved_annotations: dict[AnnotationKey, dict] = {}  # in the measured database
        self.pending_annotations: dict[AnnotationKey, dict] = {}  # to be added
        self.discarded_annotations: set[AnnotationKey] = set()  # pending ones not to add, saved ones to remove
        self.annotation_window: AnnotationManagerWindow = None

        # overview panel: a library match, a single feature or two compared features (built in a worker thread)
        self.match_canvas: FigureCanvasTkAgg = None
        self.match_toolbar: NavigationToolbar2Tk = None
        self.overview_running: bool = False
        # latest request: (request number, key, match result, feature annotations, data to reuse, settings)
        self.pending_overview: tuple = None
        # incremented for every new content of the panel, so a slow overview can't replace newer content
        self.match_panel_request: int = 0
        # (mode, feature ids, library feature id) of the overview in the panel, requesting it again doesn't replot it
        self.shown_overview: tuple = None
        self.overview_data: OverviewData = None  # data of the shown overview, reused when only settings change
        self.overview_spectra: list[PlottedSpectrum] = None  # plotted MS2 spectra, for clicks on peaks
        self.tk_hide_isotope_peaks: tk.BooleanVar = tk.BooleanVar(self.tk_root, value=False)
        self.tk_annotation_cutoff: tk.DoubleVar = tk.DoubleVar(self.tk_root, value=0.3)
        self.annotation_cutoff_job: str = None
        # user defined annotations of MS2 peaks: {("meas" | "lib", feature id): {m/z (4 decimals): annotation}}
        self.peak_annotations: dict[SpectrumKey, dict[float, str]] = {}
        # set while selection events caused by the program are pending, see _ignore_pending_table_select_events
        self.ignore_table_select: bool = False
        # annotated MS2 peaks of the measured database offered in the fragment list: {shown text: inserted value}
        self.fragment_choices: dict[str, str] = {}
        self.fragments_text_before_selection: str = ""

        # scatter plot (m/z vs RT) of the measurement
        self.scatter_canvas: FigureCanvasTkAgg = None
        self.scatter_toolbar: NavigationToolbar2Tk = None
        self.scatter_ax = None
        self.tk_show_annotations: tk.BooleanVar = tk.BooleanVar(self.tk_root, value=True)
        self.tk_rt_window: tk.DoubleVar = tk.DoubleVar(self.tk_root, value=0.2)
        # MS2 peaks below this intensity are ignored when searching fragments and losses
        self.tk_min_peak_intensity: tk.DoubleVar = tk.DoubleVar(self.tk_root, value=100)
        # zoom factor of + and - over the scatter plot
        self.scatter_zoom_factor: float = 1.5
        self.focus_before_scatter: tk.Misc = None  # gets the focus back when the mouse leaves the scatter plot
        self.tk_circle_matches: tk.BooleanVar = tk.BooleanVar(self.tk_root, value=True)
        self.highlighted_f_ids: list[int] = None  # features with the requested fragments / losses
        self.compared_f_ids: list[int] = []  # up to two features selected by shift-click
        self.last_clicked_f_id: int = None
        # pairs of features with a marked mass difference: (entry, [(feature, partner), ...]) per entry
        self.mass_difference_pairs: list[tuple[str, list[tuple[int, int]]]] = []
        self.filtered_features_marker = None  # black circles around features in the filtered table
        self.selected_feature_marker = None  # red circle around the feature selected in the table
        self.compared_features_markers: list = []  # squares around the compared features
        self.mass_difference_lines: list = []  # lines and their end caps
        self.annotated_features_markers: list = []  # green / yellow points of annotated features

        # background colors of table rows, alternating between measured features
        self.table_row_colors: list[str] = ["#E5FFCC", "#FFCCCC", "#CCFFCC", "#CCFFFF", "#CCCCFF", "#FFE5CC"]

        self._create_gui()

    # ------------------------------------------------------------------------------------------------------------------
    # GUI layout
    # ------------------------------------------------------------------------------------------------------------------
    def _create_gui(self):
        """
        Layout:
            left:   database selection, matching parameters, filters, results table, annotation buttons
            right:  overview above scatter plot (full height)
            bottom: status bar
        """
        # status bar is packed before the main area, so it keeps its space when the window shrinks
        tk_status_bar = ttk.Frame(self.tk_root, relief="sunken")
        tk_status_bar.pack(side="bottom", fill="x")
        self.tk_status_field = ttk.Label(
            tk_status_bar, text="Please load both databases to begin.", anchor="w", padding=(5, 2)
        )
        self.tk_status_field.pack(side="left", fill="x", expand=True)
        # only shown while the library is searched, see show_progress
        self.tk_progress = ttk.Progressbar(tk_status_bar, mode="determinate", maximum=100, length=250)

        tk_main_area = ttk.PanedWindow(self.tk_root, orient="horizontal")
        tk_main_area.pack(fill="both", expand=True, padx=10, pady=(10, 5))

        tk_results_pane = ttk.Frame(tk_main_area)
        tk_main_area.add(tk_results_pane, weight=1)
        self._create_database_selection(tk_results_pane)
        self._create_matching_parameters(tk_results_pane)
        self._create_filters(tk_results_pane)
        self._create_matches_table(tk_results_pane)
        self._create_annotation_buttons(tk_results_pane)

        tk_plot_area = ttk.PanedWindow(tk_main_area, orient="vertical")
        tk_main_area.add(tk_plot_area, weight=1)
        self._create_match_plot_panel(tk_plot_area)
        self._create_scatter_plot_panel(tk_plot_area)

        # the panes would otherwise be sized by the requested size of their content, so place the dividers once the
        # window is shown: results | plots at half the width, match overview | scatter plot at 55 % of the height
        def set_initial_dividers(event):
            if event.width <= 1:
                return
            tk_main_area.unbind("<Configure>")
            tk_main_area.sashpos(0, event.width // 2)
            tk_plot_area.sashpos(0, int(event.height * 0.55))

        tk_main_area.bind("<Configure>", set_initial_dividers)

    def _create_database_selection(self, parent):
        tk_database_frame = ttk.LabelFrame(parent, text="Databases", padding=(10, 5))
        tk_database_frame.pack(fill="x", pady=(0, 5))

        ttk.Label(tk_database_frame, text="Measured database").grid(row=0, column=0, sticky="w")
        ttk.Entry(tk_database_frame, textvariable=self.tk_meas_path).grid(row=0, column=1, padx=5, pady=2, sticky="ew")
        ttk.Button(tk_database_frame, text="Browse", command=self.on_browse_meas_db).grid(row=0, column=2)

        ttk.Label(tk_database_frame, text="Library database").grid(row=1, column=0, sticky="w")
        ttk.Entry(tk_database_frame, textvariable=self.tk_lib_path).grid(row=1, column=1, padx=5, pady=2, sticky="ew")
        ttk.Button(tk_database_frame, text="Browse", command=self.on_browse_lib_db).grid(row=1, column=2)
        tk_database_frame.columnconfigure(1, weight=1)

    def _create_matching_parameters(self, parent):
        tk_parameter_frame = ttk.LabelFrame(parent, text="Matching parameters", padding=(10, 5))
        tk_parameter_frame.pack(fill="x", pady=(0, 5))

        def add_entry(row: int, column: int, label: str, variable: tk.Variable):
            ttk.Label(tk_parameter_frame, text=label).grid(row=row, column=column, sticky="w", padx=(10, 4), pady=2)
            ttk.Entry(tk_parameter_frame, textvariable=variable, width=8).grid(row=row, column=column + 1, sticky="w")

        # precursor: the stricter of both tolerances is used
        add_entry(0, 0, "Precursor m/z tol. (ppm)", self.tk_max_dmz_ppm)
        add_entry(0, 2, "and (mDa)", self.tk_max_dmz_mda)
        add_entry(1, 0, "MS2 fragment m/z tol. (mDa)", self.tk_max_dmz_ms2_mda)
        add_entry(1, 2, "Min. MS2 score", self.tk_min_ms2_score)

        ttk.Label(tk_parameter_frame, text="MS2 score metric").grid(row=0, column=4, sticky="w", padx=(10, 4), pady=2)
        self.tk_ms2_score_metric: ttk.Combobox = ttk.Combobox(
            tk_parameter_frame,
            state="readonly",
            values=['cosine_fwd', 'cosine_bwd', 'cosine_sim', 'modified_cosine_greedy'],
            width=22
        )
        self.tk_ms2_score_metric.set('modified_cosine_greedy')
        self.tk_ms2_score_metric.grid(row=0, column=5, sticky="w")
        ttk.Checkbutton(tk_parameter_frame, text="Require MS2 match", variable=self.tk_require_ms2).grid(
            row=1, column=4, columnspan=2, sticky="w", padx=(10, 0)
        )

        self.tk_find_matches_button = ttk.Button(tk_parameter_frame, text="Find Matches",
                                                 command=self.on_click_find_matches)
        self.tk_find_matches_button.grid(row=0, column=6, rowspan=2, sticky="nse", padx=(20, 0))
        tk_parameter_frame.columnconfigure(6, weight=1)

    def _create_filters(self, parent):
        tk_filter_frame = ttk.LabelFrame(parent, text="Filters (applied while typing)", padding=8)
        tk_filter_frame.pack(fill="x")

        # each filter occupies a group of four grid columns: label, entry, "-", entry
        def add_label(row: int, group: int, text: str):
            ttk.Label(tk_filter_frame, text=text).grid(row=row, column=4 * group, sticky="w", padx=(10, 4), pady=2)

        def add_range_filter(row: int, group: int, label: str) -> tuple[ttk.Entry, ttk.Entry]:
            add_label(row, group, label)
            tk_min = ttk.Entry(tk_filter_frame, width=9)
            tk_min.grid(row=row, column=4 * group + 1, sticky="w")
            ttk.Label(tk_filter_frame, text="-").grid(row=row, column=4 * group + 2, padx=2)
            tk_max = ttk.Entry(tk_filter_frame, width=9)
            tk_max.grid(row=row, column=4 * group + 3, sticky="w")
            return tk_min, tk_max

        def add_filter(row: int, group: int, label: str) -> ttk.Entry:
            add_label(row, group, label)
            tk_entry = ttk.Entry(tk_filter_frame)
            tk_entry.grid(row=row, column=4 * group + 1, columnspan=3, sticky="ew")
            return tk_entry

        # left: properties of the measured features, right: properties of the matches
        self.tk_filter_min_f_id_meas, self.tk_filter_max_f_id_meas = add_range_filter(0, 0, "Meas. feature ID")
        self.tk_filter_min_rt, self.tk_filter_max_rt = add_range_filter(1, 0, "RT (min)")
        self.tk_filter_min_ccs, self.tk_filter_max_ccs = add_range_filter(2, 0, "CCS (Å²)")
        self.tk_filter_min_mz, self.tk_filter_max_mz = add_range_filter(3, 0, "m/z")
        self.tk_filter_min_intensity, self.tk_filter_max_intensity = add_range_filter(4, 0, "Σ intensity")
        self.tk_filter_min_ms2_score = add_filter(0, 1, "Min. MS2 score")
        self.tk_filter_min_n_hits_ms2 = add_filter(1, 1, "Min. #MS2 peak matches")
        self.tk_filter_subformula = add_filter(2, 1, "Formula contains")
        self.tk_filter_name = add_filter(3, 1, "Name contains (a; b)")

        tk_filter_buttons = ttk.Frame(tk_filter_frame)
        tk_filter_buttons.grid(row=4, column=4, columnspan=4, sticky="e")
        ttk.Button(tk_filter_buttons, text="m/z from compound ...", command=self.on_click_mz_from_compound).pack(
            side="left", padx=(0, 5))
        ttk.Button(tk_filter_buttons, text="Remove filters", command=self.on_click_remove_filters).pack(side="left")

        self.tk_filter_entries: list[ttk.Entry] = [
            self.tk_filter_min_f_id_meas, self.tk_filter_max_f_id_meas, self.tk_filter_name,
            self.tk_filter_min_rt, self.tk_filter_max_rt, self.tk_filter_subformula,
            self.tk_filter_min_ms2_score, self.tk_filter_min_n_hits_ms2, self.tk_filter_min_mz, self.tk_filter_max_mz,
            self.tk_filter_min_ccs, self.tk_filter_max_ccs, self.tk_filter_min_intensity, self.tk_filter_max_intensity
        ]
        for tk_entry in self.tk_filter_entries:
            tk_entry.bind("<KeyRelease>", self.on_edit_filter)

    def _create_matches_table(self, parent):
        tk_table_frame = ttk.Frame(parent)
        tk_table_frame.pack(fill="both", expand=True, pady=(5, 0))

        # column name in matches_table: (heading, width, anchor), this defines the order
        # "annotation" is not in matches_table, it shows the state of the manual annotation of the match
        columns: dict[str, tuple[str, int, str]] = {
            "feature_id_meas": ("Meas. ID", 60, "center"),
            "mz_meas": ("m/z", 80, "e"),
            "RT": ("RT (min)", 60, "e"),
            "ccs_meas": ("CCS (Å²)", 65, "e"),
            "intensity": ("Σ intensity", 80, "e"),
            "name": ("Compound name", 180, "w"),
            "formula": ("Formula", 110, "w"),
            "adduct": ("Adduct", 80, "w"),
            "ms2_score": ("MS2 score", 70, "e"),
            "n_hits_ms2": ("#MS2 peaks", 75, "e"),
            "dmz_mda": ("Δm/z (mDa)", 80, "e"),
            "dmz_ppm": ("Δm/z (ppm)", 80, "e"),
            "dccs_percent": ("ΔCCS (%)", 70, "e"),
            "annotation": ("Annotation", 85, "center"),
            "source_library": ("Source", 150, "w"),
        }
        self.table_columns: list[str] = list(columns)
        self.table_headings: dict[str, str] = {column: heading for column, (heading, _, _) in columns.items()}
        # columns of matches_table shown instead (formulas that don't fit the measured m/z are in brackets)
        self.table_display_columns: dict[str, str] = {"formula": "formula_display"}
        self.tk_matches_table: ttk.Treeview = ttk.Treeview(
            tk_table_frame,
            columns=self.table_columns,
            show="headings",
            selectmode="browse"
        )
        for column, (heading, width, anchor) in columns.items():
            self.tk_matches_table.heading(column, text=heading, command=lambda c=column: self.on_click_table_heading(c))
            self.tk_matches_table.column(column, width=width, anchor=anchor, stretch=column in ("name", "source_library"))
        self.tk_matches_table.bind("<<TreeviewSelect>>", self.on_select_table_row)
        self.tk_matches_table.bind("<ButtonRelease-1>", self.on_release_table)
        self.tk_matches_table.bind("<Control-ButtonRelease-1>", self.on_ctrl_release_table)
        self.tk_matches_table.bind("<Shift-Button-1>", self.on_shift_click_table)
        # the release of a shift+click must not show the match of the selected row (on_release_table)
        self.tk_matches_table.bind("<Shift-ButtonRelease-1>", lambda event: "break")

        # treeview cells can't be selected, so copying works with Ctrl+C (selected row) and a context menu
        self.tk_table_menu = tk.Menu(self.tk_matches_table, tearoff=False)
        self.tk_matches_table.bind("<Button-3>", self.on_right_click_table)
        self.tk_matches_table.bind("<Control-c>", lambda event: self.copy_selected_row())

        y_scroll = ttk.Scrollbar(tk_table_frame, orient="vertical", command=self.tk_matches_table.yview)
        x_scroll = ttk.Scrollbar(tk_table_frame, orient="horizontal", command=self.tk_matches_table.xview)
        self.tk_matches_table.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.tk_matches_table.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")

        tk_table_frame.rowconfigure(0, weight=1)
        tk_table_frame.columnconfigure(0, weight=1)

    def _create_annotation_buttons(self, parent):
        tk_button_frame = ttk.Frame(parent)
        tk_button_frame.pack(fill="x", pady=(5, 0))
        self.tk_manage_annotations_button = ttk.Button(
            tk_button_frame, text="Manage annotations", command=self.on_click_manage_annotations
        )
        self.tk_manage_annotations_button.pack(side="right")
        ttk.Button(tk_button_frame, text="Add library match", command=self.on_click_add_annotation).pack(
            side="right", padx=5
        )
        ttk.Button(tk_button_frame, text="Add manual annotation", command=self.on_click_add_manual_annotation).pack(
            side="right"
        )
        ttk.Button(tk_button_frame, text="Mark in-source fragment",
                   command=self.on_click_mark_in_source_fragment).pack(side="right", padx=5)

    def _create_match_plot_panel(self, parent):
        self.tk_match_plot_frame = ttk.LabelFrame(
            parent, text="Overview (select a match in the table or features in the scatter plot)"
        )
        parent.add(self.tk_match_plot_frame, weight=3)

        tk_controls = ttk.Frame(self.tk_match_plot_frame, padding=(5, 2))
        tk_controls.pack(fill="x")
        ttk.Label(tk_controls, text="Feature ID(s)").pack(side="left")
        self.tk_feature_ids_entry = ttk.Entry(tk_controls, width=14)
        self.tk_feature_ids_entry.pack(side="left", padx=(4, 0))
        self.tk_feature_ids_entry.bind("<Return>", lambda event: self.on_submit_feature_ids())
        ttk.Button(tk_controls, text="Show", command=self.on_submit_feature_ids).pack(side="left", padx=(3, 15))
        ttk.Checkbutton(tk_controls, text="Hide precursor isotope peaks (±4 Da)", variable=self.tk_hide_isotope_peaks,
                        command=self.redraw_overview).pack(side="left")
        ttk.Label(tk_controls, text="Annotation cutoff").pack(side="left", padx=(15, 3))
        ttk.Scale(tk_controls, from_=0, to=1, variable=self.tk_annotation_cutoff, length=120,
                  command=self.on_change_annotation_cutoff).pack(side="left")
        self.tk_annotation_cutoff_label = ttk.Label(tk_controls, text=f"{self.tk_annotation_cutoff.get():.2f}", width=4)
        self.tk_annotation_cutoff_label.pack(side="left", padx=(3, 0))
        ttk.Label(tk_controls, foreground="gray", text="(click an MS2 peak to annotate it)").pack(side="left", padx=10)

        self.tk_overview_plot_frame = ttk.Frame(self.tk_match_plot_frame)
        self.tk_overview_plot_frame.pack(fill="both", expand=True)

    def _create_scatter_plot_panel(self, parent):
        tk_scatter_panel = ttk.LabelFrame(
            parent, text="Scatter plot m/z vs RT (click: select in table, shift-click two features: compare MS2)"
        )
        parent.add(tk_scatter_panel, weight=2)

        tk_controls = ttk.Frame(tk_scatter_panel, padding=5)
        tk_controls.pack(fill="x")
        ttk.Checkbutton(tk_controls, text="Show annotations", variable=self.tk_show_annotations,
                        command=self.redraw_scatter_plot).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(tk_controls, text="Circle matches", variable=self.tk_circle_matches,
                        command=self.mark_filtered_features).grid(row=0, column=1, sticky="w", padx=(4, 10))
        ttk.Label(tk_controls, text="Min. peak intensity").grid(row=0, column=2, sticky="w")
        ttk.Entry(tk_controls, textvariable=self.tk_min_peak_intensity, width=10).grid(row=0, column=3, sticky="w",
                                                                                        padx=(4, 10))
        ttk.Label(tk_controls, foreground="gray", text="(fragments and losses; + / - over the plot: zoom)").grid(
            row=0, column=4, columnspan=4, sticky="w")

        # highlight features with MS2 fragments / neutral losses, entries are m/z values or formulas
        ttk.Label(tk_controls, text="Fragments (a; b)").grid(row=1, column=0, sticky="w", pady=2)
        # can be typed in, the list holds the annotated peaks of the measured database (selecting one appends it)
        self.tk_fragments_entry = ttk.Combobox(tk_controls, width=22, postcommand=self._remember_fragments_text)
        self.tk_fragments_entry.grid(row=1, column=1, sticky="w", padx=(4, 10))
        self.tk_fragments_entry.bind("<<ComboboxSelected>>", self.on_select_fragment)
        ttk.Label(tk_controls, text="Losses (a; b)").grid(row=1, column=2, sticky="w")
        self.tk_losses_entry = ttk.Entry(tk_controls, width=22)
        self.tk_losses_entry.grid(row=1, column=3, sticky="w", padx=(4, 10))
        ttk.Label(tk_controls, text="Match").grid(row=1, column=4, sticky="w")
        self.tk_highlight_mode: ttk.Combobox = ttk.Combobox(
            tk_controls, state="readonly", values=["all", "any"], width=4
        )
        self.tk_highlight_mode.set("all")
        self.tk_highlight_mode.grid(row=1, column=5, sticky="w", padx=(4, 10))
        ttk.Button(tk_controls, text="Highlight", command=self.on_click_highlight).grid(row=1, column=6, padx=3)
        ttk.Button(tk_controls, text="Remove highlight", command=self.on_click_remove_highlight).grid(row=1, column=7)
        for tk_entry in (self.tk_fragments_entry, self.tk_losses_entry):
            tk_entry.bind("<Return>", lambda event: self.on_click_highlight())

        # connect pairs of features with mass differences or adduct exchanges, matched with the precursor tolerance
        ttk.Label(tk_controls, text="Mass differences").grid(row=2, column=0, sticky="w", pady=2)
        self.tk_mass_differences_entry = ttk.Entry(tk_controls)
        self.tk_mass_differences_entry.grid(row=2, column=1, columnspan=3, sticky="ew", padx=(4, 10))
        self.tk_mass_differences_entry.bind("<Return>", lambda event: self.on_click_mark_mass_differences())
        ttk.Label(tk_controls, text="RT window (min)").grid(row=2, column=4, sticky="w")
        ttk.Entry(tk_controls, textvariable=self.tk_rt_window, width=5).grid(row=2, column=5, sticky="w", padx=(4, 10))
        ttk.Button(tk_controls, text="Mark", command=self.on_click_mark_mass_differences).grid(row=2, column=6,
                                                                                               padx=3)
        ttk.Button(tk_controls, text="Remove marks", command=self.on_click_remove_mass_differences).grid(row=2,
                                                                                                         column=7)
        ttk.Label(tk_controls, foreground="gray",
                  text="e.g. NH3; C[13] - C[12]; [M+H]+>[M+NH4]+ (adduct exchange); fragments are ions, "
                       "e.g. C5H15NO4P+").grid(
            row=3, column=0, columnspan=8, sticky="w")

        self.tk_scatter_plot_frame = ttk.Frame(tk_scatter_panel)
        self.tk_scatter_plot_frame.pack(fill="both", expand=True)

    # ------------------------------------------------------------------------------------------------------------------
    # Loading databases and closing
    # ------------------------------------------------------------------------------------------------------------------
    def on_browse_meas_db(self):
        """Connect to a measured database, reset all results and start loading its MS2 spectra."""
        if not self._confirm_discarding_unsaved_annotations():
            return
        path: str = filedialog.askopenfilename(filetypes=SQLITE_FILE_TYPES)
        if not path:
            return
        try:
            self.meas_db = FeatureManagerDB(path)
        except Exception as error:
            self.meas_db = None
            messagebox.showerror("Error", f"Failed to load measured database:\n{error}")
            return
        self.tk_meas_path.set(path)
        self.meas_db_load_number += 1

        # clear data belonging to the previous database
        self.ms2_spectra = None
        self.ms2_loading_error = None
        self.ms2_loading = True
        self.matches = {}
        self.matches_table = pd.DataFrame(columns=self.table_columns)
        self.filtered_matches_table = self.matches_table.copy()
        self._clear_matches_table()
        self.pending_annotations.clear()
        self.discarded_annotations.clear()
        self._load_saved_annotations()
        self._on_annotations_changed()

        self.match_panel_request += 1
        self.pending_overview = None
        self.shown_overview = None
        self.overview_data = self.overview_spectra = None
        self.peak_annotations = {}
        self._destroy_match_plot()

        self.highlighted_f_ids = None
        self.compared_f_ids = []
        self.last_clicked_f_id = None
        self.mass_difference_pairs = []
        self.draw_scatter_plot(keep_view=False)
        self._update_fragment_choices()

        self.tk_status_field.config(text=f"Measured database loaded, loading MS2 spectra ...: {path}")
        threading.Thread(
            target=self._ms2_spectra_worker,
            args=(self.meas_db, self.meas_db_load_number),
            daemon=True
        ).start()

    def _ms2_spectra_worker(self, meas_db: FeatureManagerDB, load_number: int):
        """Worker thread: load all MS2 spectra of the measurement."""
        try:
            spectra = meas_db.get_ms_spectra(meas_db.feature_ids, level=2)
            error_text = None
        except Exception as error:
            spectra = None
            error_text = str(error)

        try:
            self.tk_root.after(0, lambda: self._finish_ms2_loading(load_number, spectra, error_text))
        except (tk.TclError, RuntimeError):  # window was closed in the meantime
            pass

    def _finish_ms2_loading(self, load_number: int, spectra: dict | None, error_text: str | None):
        """Main thread: store the loaded MS2 spectra, unless another database was loaded in the meantime."""
        if load_number != self.meas_db_load_number:
            return
        self.ms2_loading = False
        self.ms2_spectra = spectra
        self.ms2_loading_error = error_text

        if error_text is not None:
            self.tk_status_field.config(text="Measured database loaded, but MS2 spectra could not be loaded.")
            messagebox.showerror("MS2 loading error", f"Could not load MS2 spectra:\n{error_text}")
            return
        self.tk_status_field.config(text="Measured database and MS2 spectra loaded.")

    def on_browse_lib_db(self):
        path: str = filedialog.askopenfilename(filetypes=SQLITE_FILE_TYPES)
        if not path:
            return
        try:
            self.lib_db = Library(path)
        except Exception as error:
            self.lib_db = None
            messagebox.showerror("Error", f"Failed to load library database:\n{error}")
            return
        self.tk_lib_path.set(path)
        self.tk_status_field.config(text=f"Library database loaded: {path}")
        self._clear_results_of_previous_library()
        self._offer_missing_indexes(self.lib_db)

    def _clear_results_of_previous_library(self):
        """Matches refer to library feature ids of the previous library, so they can't be shown with the new one."""
        if self.meas_db is None:
            return
        self.matches = {}
        self.matches_table = pd.DataFrame(columns=self.table_columns + ['feature_id_lib', 'formula_display',
                                                                         'default_order'])
        self.filtered_matches_table = self.matches_table.copy()
        self._clear_matches_table()
        self.mark_filtered_features()
        if self.shown_overview is not None and self.shown_overview[0] == "match":
            self.match_panel_request += 1  # drops a match overview that is still being built
            self.pending_overview = None
            self.shown_overview = None
            self.overview_data = self.overview_spectra = None
            self._destroy_match_plot()

    def _offer_missing_indexes(self, db: FeatureManagerDB):
        """
        Databases created before the indexes were added to the models load spectra much slower (full-scale library:
        search 66 s instead of 14 s). Creating them changes the file, so ask first.
        """
        try:
            missing = db.get_missing_indexes()
        except Exception:  # e.g. not an msIO database, the search will report it
            return
        if not missing:
            return
        size_gb: float = os.path.getsize(db.path_file) / 1e9
        if not messagebox.askyesno(
                "Speed up the library",
                f"The library has no index for loading spectra, which makes the search several times slower.\n\n"
                f"Add the index now? This is done once and changes the file: it takes a while for large libraries "
                f"(about 20 s per 5 GB) and makes the file about 20 % larger (currently {size_gb:.1f} GB)."
        ):
            return
        self.tk_find_matches_button.config(state="disabled")
        self.show_progress("Creating the index of the library ...", indeterminate=True)

        def worker():
            try:
                db.create_missing_indexes()
                error_text = None
            except Exception as error:
                error_text = f"{type(error).__name__}: {error}"
            self._call_in_main_thread(lambda: finish(error_text))

        def finish(error_text: str | None):
            self.tk_find_matches_button.config(state="normal")
            self.hide_progress()
            if error_text is not None:
                messagebox.showerror("Index", f"Could not create the index:\n{error_text}")
                self.tk_status_field.config(text="Library loaded without index.")
            else:
                self.tk_status_field.config(text="Index of the library created.")

        threading.Thread(target=worker, daemon=True).start()

    def on_close(self):
        """Closing the main window: warn about annotation changes that were not written to the database."""
        if self._confirm_discarding_unsaved_annotations():
            self.tk_root.destroy()

    # ------------------------------------------------------------------------------------------------------------------
    # Matching
    # ------------------------------------------------------------------------------------------------------------------
    def on_click_find_matches(self):
        """
        Validate the matching parameters and search the library for all measured features in a worker thread
        (progress in the status bar), the table is filled when the search is done.
        """
        if self.meas_db is None or self.lib_db is None:
            messagebox.showerror("Error", "Please select both databases.")
            return

        require_ms2: bool = self.tk_require_ms2.get()
        # without spectra, all MS2 scores would be missing
        if self.ms2_loading:
            messagebox.showinfo("MS2 spectra", "MS2 spectra are still loading. Try again shortly.")
            return
        if require_ms2 and (self.ms2_loading_error or self.ms2_spectra is None):
            messagebox.showerror("MS2 spectra", "MS2 spectra could not be loaded.")
            return

        try:
            max_dmz_ppm: float = self.tk_max_dmz_ppm.get()
            max_dmz_da: float = self.tk_max_dmz_mda.get() * 1e-3
            max_ms2_dmz_da: float = self.tk_max_dmz_ms2_mda.get() * 1e-3
            min_ms2_score: float = self.tk_min_ms2_score.get()
        except tk.TclError:
            messagebox.showerror("Invalid input", "Matching parameters must be numbers.")
            return
        metric: str = self.tk_ms2_score_metric.get()

        if max_dmz_ppm <= 0 or max_dmz_da <= 0 or max_ms2_dmz_da <= 0:
            messagebox.showerror("Invalid input", "Tolerances must be greater than 0.")
            return
        if not 0 <= min_ms2_score <= 1:
            messagebox.showerror("Invalid input", "Please enter a MS2 score between 0 and 1.")
            return

        self.tk_find_matches_button.config(state="disabled")
        self.show_progress("Searching the library ...")
        search_kwargs: dict = dict(
            mzs=self.meas_db.mzs,
            max_dmz_ppm=max_dmz_ppm,
            max_dmz_da=max_dmz_da,
            max_ms2_dmz_da=max_ms2_dmz_da,
            metric=metric,
            min_ms2_score=min_ms2_score,
            return_nhits_ms2=True,
            # empty dict: scores are missing, but find_matches still returns the matches as dicts
            ms2_spectra=self.ms2_spectra or {},
            require_ms2=require_ms2,
        )
        threading.Thread(
            target=self._find_matches_worker,
            args=(self.meas_db, self.lib_db, search_kwargs, self.meas_db_load_number),
            daemon=True
        ).start()

    def _find_matches_worker(self, meas_db: FeatureManagerDB, lib_db: Library, search_kwargs: dict, load_number: int):
        """Worker thread: search the library and build the results table."""
        last_percent: list[int] = [-1]

        def report_progress(n_done: int, n_total: int):
            # only hand over changes of whole percents, the main thread doesn't need every feature
            percent = int(n_done / max(n_total, 1) * 100)
            if percent != last_percent[0]:
                last_percent[0] = percent
                self._call_in_main_thread(lambda: self.update_progress(percent, f"Scoring MS2 matches ... {percent} %"))

        matches, table, error_text = None, None, None
        try:
            matches = lib_db.find_matches(**search_kwargs, progress_callback=report_progress)
            self._call_in_main_thread(lambda: self.update_progress(100, "Building the results table ..."))
            table = self._matches_to_table(matches, meas_db, lib_db, search_kwargs['max_dmz_ppm'],
                                           search_kwargs['max_dmz_da'])
        except Exception as error:
            error_text = f"{type(error).__name__}: {error}"
        self._call_in_main_thread(lambda: self._finish_find_matches(load_number, lib_db, matches, table, error_text))

    def _call_in_main_thread(self, func: Callable[[], None]):
        """Hand a function from a worker thread to the tk main loop."""
        try:
            self.tk_root.after(0, func)
        except (tk.TclError, RuntimeError):  # window was closed in the meantime
            pass

    def _finish_find_matches(self, load_number: int, lib_db: Library, matches: dict | None, table: pd.DataFrame | None,
                             error_text: str | None):
        """Main thread: show the results, unless another database was loaded during the search."""
        self.tk_find_matches_button.config(state="normal")
        self.hide_progress()
        if load_number != self.meas_db_load_number or lib_db is not self.lib_db:
            self.tk_status_field.config(text="Search results dropped, a different database was loaded meanwhile.")
            return
        if error_text is not None:
            messagebox.showerror("Error", f"Search failed:\n{error_text}")
            self.tk_status_field.config(text="Failed while looking for matches.")
            return

        self.matches, self.matches_table = matches, table
        self.shown_overview = None  # scores may have changed
        self.apply_filters(show_empty_message=False)
        count: int = len(self.matches_table)
        self.tk_status_field.config(text=f"Found {count} matches.")
        if count == 0:
            messagebox.showinfo("Search complete", "No matches found.")

    def show_progress(self, text: str, indeterminate: bool = False):
        """Show the progress bar in the status bar, indeterminate: moving bar for work without known progress."""
        self.tk_status_field.config(text=text)
        self.tk_progress.config(value=0, mode="indeterminate" if indeterminate else "determinate")
        self.tk_progress.pack(side="right", padx=5, before=self.tk_status_field)
        if indeterminate:
            self.tk_progress.start(20)

    def update_progress(self, percent: float, text: str):
        self.tk_progress.config(value=percent)
        self.tk_status_field.config(text=text)

    def hide_progress(self):
        self.tk_progress.stop()
        self.tk_progress.pack_forget()

    def _matches_to_table(
            self,
            matches: dict[int, list[dict]],
            meas_db: FeatureManagerDB,
            lib_db: Library,
            max_dmz_ppm: float,
            max_dmz_da: float
    ) -> pd.DataFrame:
        """
        Turn the output of Library.find_matches into the table shown in the GUI.

        Adds properties of the measured features (m/z, RT in minutes, CCS, summed intensity), the CCS difference
        to the library entry and the formula to display (in brackets if it doesn't fit the measured m/z).
        Runs in the search worker, so the databases are passed in instead of using the current ones.
        """
        if not any(matches.values()):  # lib_match_result_to_table can't handle that
            return pd.DataFrame(columns=self.table_columns + ['feature_id_lib', 'formula_display', 'default_order'])
        table: pd.DataFrame = lib_match_result_to_table(matches)
        f_ids_meas: pd.Series = table.feature_id_meas
        table.loc[:, 'ms2_score'] = table['ms2_score'].round(3)
        table.loc[:, 'dmz_mda'] = table['dmz_mda'].round(1).abs()
        table.loc[:, 'dmz_ppm'] = table['dmz_ppm'].round(1).abs()
        table['RT'] = (f_ids_meas.map(meas_db.retention_times_in_seconds) / 60).round(2)
        table['intensity'] = pd.to_numeric(f_ids_meas.map(meas_db.summed_intensities), errors='coerce').round(2)

        # CCS of 0 means not measured
        ccs_meas = pd.to_numeric(f_ids_meas.map(meas_db.collisional_cross_sections), errors='coerce')
        ccs_lib = pd.to_numeric(table.feature_id_lib.map(lib_db.collisional_cross_sections), errors='coerce')
        table['ccs_meas'] = ccs_meas.where(ccs_meas > 0).round(1)
        table['dccs_percent'] = ((ccs_meas - ccs_lib) / ccs_lib * 100).where((ccs_meas > 0) & (ccs_lib > 0)).round(2)

        mzs_meas = f_ids_meas.map(meas_db.mzs)
        has_formula = ~table.formula.isin(['', 'None', 'nan'])
        fits_mz = self._formula_fits_mz(table.formula, table.adduct, mzs_meas.to_numpy(), max_dmz_ppm, max_dmz_da)
        table['formula_display'] = np.where(
            has_formula, np.where(fits_mz, table.formula, '(' + table.formula + ')'), ''
        )
        table['mz_meas'] = mzs_meas.round(4)
        # order of the search result, restored when the sorting is removed (see _sort_matches)
        table['default_order'] = np.arange(len(table))
        return table

    @staticmethod
    def _formula_fits_mz(
            formulas: pd.Series,
            adducts: pd.Series,
            mzs_meas: np.ndarray,
            max_dmz_ppm: float,
            max_dmz_da: float
    ) -> np.ndarray:
        """
        Whether the m/z calculated from formula and adduct matches the measured m/z within the precursor tolerance
        (stricter of both, as in the search). Formulas or adducts that can't be parsed don't fit.
        """
        def calculate_mz(formula: str, adduct: str) -> float:
            try:
                return Adduct(adduct).mass_to_mz(CompoundDict(formula).mass)
            except Exception:
                return np.nan

        # many matches share formula and adduct
        mzs_calculated: dict[tuple[str, str], float] = {
            key: calculate_mz(*key) for key in set(zip(formulas, adducts))
        }
        mzs_calc = np.array([mzs_calculated[key] for key in zip(formulas, adducts)], dtype=float)
        tolerances = np.minimum(max_dmz_da, max_dmz_ppm * 1e-6 * mzs_meas)
        # nan (not parsable) compares as False
        return np.abs(mzs_calc - mzs_meas) <= tolerances

    # ------------------------------------------------------------------------------------------------------------------
    # Filters
    # ------------------------------------------------------------------------------------------------------------------
    def apply_filters(self, show_empty_message: bool = True, show_dialogs: bool = True):
        """
        Filter matches_table by all filter fields and show the result in the table.

        Invalid filter values leave the table unchanged. They are reported in a dialog,
        or in the status bar if show_dialogs is False (used while typing).
        """
        def report_invalid(title: str, text: str):
            # while typing, half-entered values are expected, so don't open a dialog for every key stroke
            if show_dialogs:
                messagebox.showerror(title, text)
            else:
                self.tk_status_field.config(text=f"Filters not updated: {text}")

        if self.matches_table.empty:
            self._set_filtered_matches_table(self.matches_table.copy())
            return

        table: pd.DataFrame = self.matches_table
        mask = np.ones(table.shape[0], dtype=bool)

        # numeric filters: (entry, column to compare, type, is minimum, error message)
        # features without a value (e.g. no CCS) are removed by filters on that column
        numeric_filters = [
            (self.tk_filter_min_f_id_meas, table.feature_id_meas, int, True, 'Minimum ID must be a whole number.'),
            (self.tk_filter_max_f_id_meas, table.feature_id_meas, int, False, 'Maximum ID must be a whole number.'),
            (self.tk_filter_min_rt, table.RT, float, True, 'Minimum RT must be a number.'),
            (self.tk_filter_max_rt, table.RT, float, False, 'Maximum RT must be a number.'),
            (self.tk_filter_min_mz, table.mz_meas, float, True, 'Minimum m/z must be a number.'),
            (self.tk_filter_max_mz, table.mz_meas, float, False, 'Maximum m/z must be a number.'),
            (self.tk_filter_min_ms2_score, table.ms2_score, float, True, 'Minimum MS2 score must be a number.'),
            (self.tk_filter_min_n_hits_ms2, table.n_hits_ms2, int, True,
             'Minimum #MS2 peak matches must be a whole number.'),
            (self.tk_filter_min_ccs, table.ccs_meas, float, True, 'Minimum CCS must be a number.'),
            (self.tk_filter_max_ccs, table.ccs_meas, float, False, 'Maximum CCS must be a number.'),
            (self.tk_filter_min_intensity, table.intensity, float, True, 'Minimum intensity must be a number.'),
            (self.tk_filter_max_intensity, table.intensity, float, False, 'Maximum intensity must be a number.'),
        ]
        for tk_entry, values, value_type, is_minimum, error_message in numeric_filters:
            if not (text := tk_entry.get().strip()):
                continue
            try:
                limit = value_type(text)
            except ValueError:
                report_invalid('Value Error', error_message)
                return
            mask &= (values >= limit) if is_minimum else (values <= limit)

        # several name substrings can be separated by ";", all of them must be contained
        if name_substrings := self.tk_filter_name.get().strip():
            for substring in name_substrings.split(';'):
                mask &= table.name.astype(str).str.contains(substring.strip(), case=False, na=False, regex=False)

        # formula must contain at least as many atoms of each element as the given subformula
        if subformula := self.tk_filter_subformula.get().strip():
            try:
                composition_sub: dict[str, int] = CompoundDict(subformula).composition
                mask &= np.array([
                    all(CompoundDict(formula).composition.get(element, 0) >= count
                        for element, count in composition_sub.items())
                    for formula in table.formula
                ])
            except Exception as error:
                report_invalid('Invalid Formula', f"Invalid formula: {error}")
                return

        self._set_filtered_matches_table(table.loc[mask, :])
        if show_empty_message and self.filtered_matches_table.empty:
            self.tk_status_field.config(text="No matches meet the current filters.")
        elif not self.filtered_matches_table.empty:
            self.tk_status_field.config(
                text=f"Showing {len(self.filtered_matches_table)} of {len(self.matches_table)} matches.")

    def _get_filter_values(self) -> tuple[str, ...]:
        return tuple(tk_entry.get().strip() for tk_entry in self.tk_filter_entries)

    def on_edit_filter(self, event=None):
        """Key release in a filter field: apply the filters once typing pauses."""
        # keys like Tab or arrows don't change the filters
        if self._get_filter_values() == self.last_filter_values:
            return
        # restart the timer on every key stroke, so filtering only happens once typing pauses
        if self.live_filter_job is not None:
            self.tk_root.after_cancel(self.live_filter_job)
        self.live_filter_job = self.tk_root.after(self.live_filter_delay_ms, self._apply_live_filters)

    def _apply_live_filters(self):
        self.live_filter_job = None
        self.last_filter_values = self._get_filter_values()
        self.apply_filters(show_dialogs=False)

    def get_precursor_tolerance_da(self, mz: float) -> float:
        """Precursor m/z tolerance in Da at mz: the stricter of the ppm and mDa tolerance (as in the search)."""
        return min(self.tk_max_dmz_mda.get() * 1e-3, self.tk_max_dmz_ppm.get() * 1e-6 * mz)

    def set_mz_filter(self, mz: float):
        """Filter the results to measured features at mz (± precursor tolerance)."""
        tolerance: float = self.get_precursor_tolerance_da(mz)
        for tk_entry, value in ((self.tk_filter_min_mz, mz - tolerance), (self.tk_filter_max_mz, mz + tolerance)):
            tk_entry.delete(0, tk.END)
            tk_entry.insert(0, f"{value:.4f}")
        self.last_filter_values = self._get_filter_values()
        self.apply_filters()

    def on_click_mz_from_compound(self):
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        MzFromCompoundDialog(self)

    def on_click_remove_filters(self):
        """Empty all filter fields and show all matches."""
        for tk_entry in self.tk_filter_entries:
            tk_entry.delete(0, tk.END)
        self.last_filter_values = self._get_filter_values()
        self._set_filtered_matches_table(self.matches_table)
        self.tk_status_field.config(text=f"Showing all {len(self.matches_table)} matches.")

    # ------------------------------------------------------------------------------------------------------------------
    # Results table
    # ------------------------------------------------------------------------------------------------------------------
    def _clear_matches_table(self):
        self.tk_matches_table.delete(*self.tk_matches_table.get_children())

    @staticmethod
    def _library_annotation_key(f_id_meas: int, f_id_lib: int) -> AnnotationKey:
        return "library", int(f_id_meas), str(int(f_id_lib))

    def _set_filtered_matches_table(self, table: pd.DataFrame):
        """Sort (by the clicked columns) and show the filtered matches, the selected match stays selected."""
        selected = self._get_selected_match()
        if not table.empty:
            table = self._sort_matches(table, self.sort_keys)
        self.filtered_matches_table = table.reset_index(drop=True)
        self._show_matches_in_table(self.filtered_matches_table)
        self.mark_filtered_features()

        if selected is not None:
            is_selected = ((self.filtered_matches_table.feature_id_meas == selected.feature_id_meas)
                           & (self.filtered_matches_table.feature_id_lib == selected.feature_id_lib))
            if is_selected.any():
                # restoring the selection must not change the overview (e.g. a feature shown from the scatter plot)
                self._ignore_pending_table_select_events()
                item: str = self.tk_matches_table.get_children()[int(np.argmax(is_selected))]
                self.tk_matches_table.selection_set(item)
                self.tk_matches_table.see(item)

    def _sort_matches(self, table: pd.DataFrame, sort_keys: list[tuple[str, bool]]) -> pd.DataFrame:
        """
        Sort matches by (column, descending) pairs, the first one is the primary key. Text is sorted case-insensitive,
        missing values last, ties keep the default order (as found by the search).
        """
        table = table.sort_values("default_order", kind="mergesort")
        # stable sorts from the last key to the first one give the sorting by several keys
        for column, descending in reversed(sort_keys):
            if column == "annotation":  # not a column of the matches table
                values = pd.Series([self._get_annotation_state(self._library_annotation_key(f_id_meas, f_id_lib))
                                    for f_id_meas, f_id_lib in zip(table.feature_id_meas, table.feature_id_lib)],
                                   index=table.index)
            else:
                values = table[column]
            if values.dtype == object:
                values = values.astype(str).str.lower()
            table = table.loc[values.sort_values(ascending=not descending, na_position="last", kind="mergesort").index]
        return table

    def on_click_table_heading(self, column: str, secondary: bool = False):
        """
        Click on a column header: sort by it, ascending -> descending -> default order.
        Ctrl+click (secondary): add the column as further sort key (same cycle for that key).
        """
        keys: list[tuple[str, bool]] = self.sort_keys
        position = next((i for i, (c, _) in enumerate(keys) if c == column), None)
        if not secondary or not keys:
            if position == 0:  # primary key: reverse, then back to the default order
                self.sort_keys = [(column, True)] + keys[1:] if not keys[0][1] else []
            else:
                self.sort_keys = [(column, False)]
        elif position is None:
            self.sort_keys = keys + [(column, False)]
        elif not keys[position][1]:
            self.sort_keys = keys[:position] + [(column, True)] + keys[position + 1:]
        else:
            self.sort_keys = keys[:position] + keys[position + 1:]

        numbered: bool = len(self.sort_keys) > 1
        arrows = {c: (" ▼" if descending else " ▲") + (str(i + 1) if numbered else "")
                  for i, (c, descending) in enumerate(self.sort_keys)}
        for c, heading in self.table_headings.items():
            self.tk_matches_table.heading(c, text=heading + arrows.get(c, ""))
        self._set_filtered_matches_table(self.filtered_matches_table)

    def on_ctrl_release_table(self, event):
        """Ctrl+click on a column header: further sort key."""
        if self.tk_matches_table.identify_region(event.x, event.y) != "heading":
            return None
        column: str = self.table_columns[int(self.tk_matches_table.identify_column(event.x)[1:]) - 1]
        self.on_click_table_heading(column, secondary=True)
        return "break"  # don't invoke the normal sorting of the header

    def on_shift_click_table(self, event):
        """Shift+click on a row: add its measured feature to the comparison (like shift+click in the scatter plot)."""
        if not (item := self.tk_matches_table.identify_row(event.y)):
            return None
        row = self.filtered_matches_table.iloc[self.tk_matches_table.index(item)]
        self.add_feature_to_comparison(int(row.feature_id_meas))
        return "break"  # keep the selection

    def _show_matches_in_table(self, table: pd.DataFrame):
        """Replace the rows of the results table, rows of the same measured feature share a background color."""
        self._clear_matches_table()

        f_ids_meas = sorted(table['feature_id_meas'].unique())
        n_colors: int = len(self.table_row_colors)
        for color_idx, f_id_meas in enumerate(f_ids_meas):
            self.tk_matches_table.tag_configure(
                f"feature_{f_id_meas}", background=self.table_row_colors[color_idx % n_colors]
            )

        # values of all rows at once (iterrows takes seconds for tens of thousands of matches)
        source_columns: list[str] = [self.table_display_columns.get(column, column) for column in self.table_columns]
        cells = pd.DataFrame({
            column: table[column] if column in table.columns else pd.Series("", index=table.index)
            for column in source_columns
        }).astype(object)
        rows: list[list] = cells.where(cells.notna(), "").to_numpy().tolist()
        idx_annotation: int = self.table_columns.index("annotation")
        for values, f_id_meas, f_id_lib in zip(rows, table.feature_id_meas, table.feature_id_lib):
            values[idx_annotation] = self._get_annotation_state(self._library_annotation_key(f_id_meas, f_id_lib))
            self.tk_matches_table.insert("", 'end', values=values, tags=(f"feature_{f_id_meas}",))

    def _copy_to_clipboard(self, text: str, description: str):
        self.tk_root.clipboard_clear()
        self.tk_root.clipboard_append(text)
        self.tk_status_field.config(text=f"Copied {description} to the clipboard.")

    def _get_row_text(self, item: str) -> str:
        return "\t".join(str(value) for value in self.tk_matches_table.item(item, "values"))

    def copy_selected_row(self):
        if selection := self.tk_matches_table.selection():
            self._copy_to_clipboard(self._get_row_text(selection[0]), "the selected row")

    def copy_table(self):
        """Copy the shown (filtered and sorted) results table as tab separated text (e.g. to paste into Excel)."""
        lines = ["\t".join(self.table_headings[column] for column in self.table_columns)]
        lines += [self._get_row_text(item) for item in self.tk_matches_table.get_children()]
        self._copy_to_clipboard("\n".join(lines), f"{len(lines) - 1} rows")

    def on_right_click_table(self, event):
        """Select the row under the mouse and offer to copy the cell, the row or the whole table."""
        item: str = self.tk_matches_table.identify_row(event.y)
        if not item:
            return
        self.tk_matches_table.selection_set(item)
        column: str = self.table_columns[int(self.tk_matches_table.identify_column(event.x)[1:]) - 1]
        value: str = str(self.tk_matches_table.set(item, column))
        self.tk_table_menu.delete(0, "end")
        self.tk_table_menu.add_command(
            label=f"Copy {self.table_headings[column]}: {value if len(value) <= 40 else value[:37] + '...'}",
            command=lambda: self._copy_to_clipboard(value, f'"{value}"')
        )
        self.tk_table_menu.add_command(label="Copy row", command=self.copy_selected_row)
        self.tk_table_menu.add_command(label="Copy table", command=self.copy_table)
        self.tk_table_menu.tk_popup(event.x_root, event.y_root)

    def _get_selected_match(self) -> pd.Series | None:
        """Row of filtered_matches_table that is selected in the results table."""
        selection: tuple[str, ...] = self.tk_matches_table.selection()
        if not selection:
            return None
        try:
            return self.filtered_matches_table.iloc[self.tk_matches_table.index(selection[0])]
        except IndexError:
            return None

    def on_select_table_row(self, event=None):
        """Mark the selected feature in the scatter plot and plot the match overview."""
        # also fired with an empty selection when the table is refilled
        if (match := self._get_selected_match()) is None:
            return
        self.mark_selected_feature(int(match.feature_id_meas))
        if self.ignore_table_select:  # selected by the program (scatter click, restored selection): keep the overview
            return
        # e.g. selected again after sorting or filtering
        if ("match", (int(match.feature_id_meas),), int(match.feature_id_lib)) != self.shown_overview:
            self.clear_comparison()
            self.request_overview("match", [int(match.feature_id_meas)], int(match.feature_id_lib))

    def _ignore_pending_table_select_events(self):
        """
        Selection events of the table caused by the program (e.g. deleting and restoring the selected row) arrive
        later through the event queue, possibly several of them: ignore all of them until the queue is processed.
        """
        self.ignore_table_select = True
        self.tk_root.after_idle(lambda: setattr(self, "ignore_table_select", False))

    def on_release_table(self, event):
        """Click on the already selected row (no selection event): show its match again, e.g. after a feature."""
        if (item := self.tk_matches_table.identify_row(event.y)) and item in self.tk_matches_table.selection():
            self.on_select_table_row()

    # ------------------------------------------------------------------------------------------------------------------
    # Manual annotations
    # ------------------------------------------------------------------------------------------------------------------
    def _load_saved_annotations(self):
        """Read the manual annotations (all kinds) stored in the measured database."""
        self.saved_annotations = {}
        for candidate in self.meas_db.get_manual_annotations():
            kind: str = FeatureManagerDB.get_manual_annotation_kind(candidate)
            reference: int | None = FeatureManagerDB.get_referenced_feature_id(candidate)
            ref: str = f"{candidate.name_sirius}|{candidate.adduct_sirius}" if kind == "manual" else str(reference)
            key: AnnotationKey = (kind, candidate.feature_id, ref)
            if key in self.saved_annotations:  # stored more than once
                self.saved_annotations[key]["candidate_ids"].append(candidate.id)
                continue
            self.saved_annotations[key] = dict(
                kind=kind,
                name=candidate.name_sirius,
                formula=candidate.formula_sirius,
                adduct=candidate.adduct_sirius,
                smiles=candidate.smiles,
                ms2_score=candidate.confidence_score,
                reference=reference,
                candidate_ids=[candidate.id],
            )

    def _get_annotation_state(self, key: AnnotationKey) -> str:
        """State of a manual annotation: "in database", "to remove", "to add" or "" (none)."""
        if key in self.saved_annotations:
            return "to remove" if key in self.discarded_annotations else "in database"
        if key in self.pending_annotations and key not in self.discarded_annotations:
            return "to add"
        return ""

    def get_unsaved_annotation_changes(self) -> tuple[list[AnnotationKey], list[AnnotationKey]]:
        """Annotations to add to and to remove from the measured database."""
        to_add = [key for key in self.pending_annotations if key not in self.discarded_annotations]
        to_remove = [key for key in self.saved_annotations if key in self.discarded_annotations]
        return to_add, to_remove

    def _confirm_discarding_unsaved_annotations(self) -> bool:
        to_add, to_remove = self.get_unsaved_annotation_changes()
        if not (to_add or to_remove):
            return True
        return messagebox.askyesno(
            "Unsaved annotations",
            f"{len(to_add)} annotation(s) to add and {len(to_remove)} to remove have not been written to the "
            f"database. Discard these changes?"
        )

    def add_pending_annotation(self, key: AnnotationKey, info: dict) -> bool:
        """
        Add an annotation to the local list (not yet written to the database).

        Returns False if the annotation is already in the database (it is kept, if it was marked for removal).
        """
        self.discarded_annotations.discard(key)
        if key in self.saved_annotations:
            self._on_annotations_changed()
            return False
        self.pending_annotations[key] = info
        self._on_annotations_changed()
        return True

    def _get_current_feature(self) -> int | None:
        """Feature to prefill annotation dialogs: selected in the table, else last clicked in the scatter plot."""
        if (match := self._get_selected_match()) is not None:
            return int(match.feature_id_meas)
        if self.last_clicked_f_id is not None:
            return self.last_clicked_f_id
        return self.compared_f_ids[-1] if self.compared_f_ids else None

    def on_click_add_annotation(self):
        """Add the selected library match to the local list of annotations."""
        if (match := self._get_selected_match()) is None:
            messagebox.showerror("Error", "Please select a match row first.")
            return
        f_id_meas, f_id_lib = int(match.feature_id_meas), int(match.feature_id_lib)
        match_result: dict = next(m for m in self.matches[f_id_meas] if m['feature_id'] == f_id_lib)
        added: bool = self.add_pending_annotation(self._library_annotation_key(f_id_meas, f_id_lib), dict(
            kind="library",
            name=match_result.get('name'),
            formula=match_result.get('formula'),
            adduct=match_result.get('adduct'),
            ms2_score=match_result.get('ms2_score'),
            reference=f_id_lib,
            library=self.lib_db,
        ))
        if added:
            self.tk_status_field.config(
                text=f"Added library feature {f_id_lib} for feature {f_id_meas} to the annotations. "
                     f"Use \"Manage annotations\" to write them to the database."
            )
        else:
            self.tk_status_field.config(text=f"Library feature {f_id_lib} is already annotated for feature {f_id_meas}.")

    def on_click_add_manual_annotation(self):
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        ManualAnnotationDialog(self, self._get_current_feature())

    def on_click_mark_in_source_fragment(self):
        """Open the in-source fragment dialog, prefilled with the compared features (lower m/z = fragment)."""
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        if len(self.compared_f_ids) == 2:
            f_id_fragment, f_id_parent = sorted(self.compared_f_ids, key=lambda f_id: self.meas_db.mzs[f_id])
        else:
            f_id_fragment, f_id_parent = self._get_current_feature(), None
        InSourceFragmentDialog(self, f_id_fragment, f_id_parent)

    def toggle_annotations(self, keys: list[AnnotationKey], keep: bool = None):
        """Mark annotations to keep or discard (keep=None toggles each of them)."""
        for key in keys:
            discard = (key not in self.discarded_annotations) if keep is None else not keep
            if discard:
                self.discarded_annotations.add(key)
            else:
                self.discarded_annotations.discard(key)
        self._on_annotations_changed()

    def _write_annotation(self, key: AnnotationKey):
        """Add a pending annotation to the measured database."""
        kind, f_id_meas, _ = key
        info: dict = self.pending_annotations[key]
        if kind == "library":
            self.meas_db.add_annotations_from_library(
                library=info['library'], f_id_meas=f_id_meas, f_id_lib=info['reference'], ms2_score=info['ms2_score']
            )
        elif kind == "manual":
            self.meas_db.add_manual_annotation(
                f_id_meas, name=info['name'], formula=info['formula'], adduct=info['adduct'], smiles=info['smiles']
            )
        elif kind == "in_source_fragment":
            self.meas_db.add_in_source_fragment_annotation(f_id_meas, info['reference'], note=info['note'])
        else:
            raise ValueError(f'unknown annotation kind {kind}')

    def write_annotations_to_db(self, parent: tk.Misc = None):
        """Add the new and remove the discarded annotations in the measured database (after confirmation)."""
        to_add, to_remove = self.get_unsaved_annotation_changes()
        if not (to_add or to_remove):
            messagebox.showinfo("Annotations", "There are no changes to write.", parent=parent)
            return
        if not messagebox.askyesno(
                "Write annotations",
                f"Add {len(to_add)} and remove {len(to_remove)} annotation(s) in the measured database?",
                parent=parent
        ):
            return

        errors: list[str] = []
        failed_additions: set[AnnotationKey] = set()
        n_removed: int = 0
        for key in to_add:
            try:
                self._write_annotation(key)
            except Exception as error:
                failed_additions.add(key)
                errors.append(f"adding {ANNOTATION_KIND_LABELS[key[0]]} to feature {key[1]}: {error}")
        for key in to_remove:
            try:
                self.meas_db.remove_annotations_by_id(self.saved_annotations[key]['candidate_ids'])
                n_removed += 1
            except Exception as error:
                errors.append(f"removing {ANNOTATION_KIND_LABELS[key[0]]} from feature {key[1]}: {error}")

        # failed additions stay in the local list, everything else is now reflected by the database
        self.pending_annotations = {
            key: value for key, value in self.pending_annotations.items() if key in failed_additions
        }
        self.discarded_annotations.clear()
        self._load_saved_annotations()
        self._on_annotations_changed()
        if self.tk_show_annotations.get():
            self.redraw_scatter_plot()

        if errors:
            messagebox.showerror("Annotation Error", "\n".join(errors), parent=parent)
        self.tk_status_field.config(
            text=f"Wrote annotations to the database: {len(to_add) - len(failed_additions)} added, {n_removed} removed."
        )

    def _on_annotations_changed(self):
        """Update everything that shows annotation states: table column, button text and the manager window."""
        for item, row in zip(self.tk_matches_table.get_children(), self.filtered_matches_table.itertuples()):
            self.tk_matches_table.set(item, "annotation", self._get_annotation_state(
                self._library_annotation_key(row.feature_id_meas, row.feature_id_lib)
            ))
        n_changes: int = sum(len(keys) for keys in self.get_unsaved_annotation_changes())
        self.tk_manage_annotations_button.config(
            text=f"Manage annotations ({n_changes} unsaved)" if n_changes else "Manage annotations"
        )
        if self.annotation_window is not None:
            self.annotation_window.refresh()
        self.mark_annotated_features()

    def on_click_manage_annotations(self):
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        if self.annotation_window is None:
            self.annotation_window = AnnotationManagerWindow(self)
        else:
            self.annotation_window.lift()

    # ------------------------------------------------------------------------------------------------------------------
    # Match overview panel
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def _pack_canvas(canvas: FigureCanvasTkAgg):
        """Fill the parent with the canvas, the figure is scaled to the available space."""
        tk_widget = canvas.get_tk_widget()
        # without this, the canvas requests the full figure size and pushes the dividers of the panes
        tk_widget.config(width=1, height=1)
        tk_widget.pack(fill=tk.BOTH, expand=True)

    def _show_in_match_panel(self, fig: Figure, title: str):
        """Replace the content of the overview panel."""
        self._destroy_match_plot()
        self.tk_match_plot_frame.config(text=title)
        self.match_canvas = FigureCanvasTkAgg(fig, master=self.tk_overview_plot_frame)
        self.match_toolbar = NavigationToolbar2Tk(self.match_canvas, self.tk_overview_plot_frame)
        self.match_toolbar.update()
        self.match_toolbar.pack()
        self._pack_canvas(self.match_canvas)
        self.match_canvas.mpl_connect('button_press_event', self.on_click_overview)
        self.match_canvas.draw()

    def _destroy_match_plot(self):
        if self.match_canvas is not None:
            self.match_canvas.get_tk_widget().destroy()
            self.match_canvas = None
        if self.match_toolbar is not None:
            self.match_toolbar.destroy()
            self.match_toolbar = None

    def _get_overview_settings(self) -> OverviewSettings | None:
        try:
            ms2_tolerance_da: float = self.tk_max_dmz_ms2_mda.get() * 1e-3
            ms1_tolerance_da: float = self.tk_max_dmz_mda.get() * 1e-3
            annotation_cutoff: float = self.tk_annotation_cutoff.get()
        except tk.TclError:
            messagebox.showerror("Invalid Tolerance", "Check the m/z tolerances.")
            return None
        return OverviewSettings(
            hide_isotope_peaks=self.tk_hide_isotope_peaks.get(),
            annotation_cutoff=annotation_cutoff,
            ms2_tolerance_da=ms2_tolerance_da,
            ms1_tolerance_da=ms1_tolerance_da,
            # copy, the worker reads it while the user may edit annotations
            peak_annotations={key: dict(annotations) for key, annotations in self.peak_annotations.items()},
        )

    def _get_feature_annotation(self, f_id: int) -> dict | None:
        """Manual annotation of a measured feature (pending or saved): manual > library match > in-source fragment."""
        annotations: dict[str, dict] = {}
        for collection in (self.saved_annotations, self.pending_annotations):
            for key, info in collection.items():
                if key[1] == f_id and key not in self.discarded_annotations:
                    annotations.setdefault(key[0], info)
        for kind in ("manual", "library", "in_source_fragment"):
            if kind in annotations:
                return annotations[kind]
        return None

    def request_overview(self, mode: str, f_ids: list[int], f_id_lib: int = None, reuse_data: bool = False):
        """
        Show an overview in the panel: "match" (feature f_ids[0] and library feature f_id_lib), "feature" (f_ids[0])
        or "comparison" (both f_ids).

        Gathering the data and building the figure are slow, so both run in a worker thread. Only one worker runs at
        a time; if another overview is requested meanwhile, only the latest one is shown afterwards.
        With reuse_data, only the settings changed (cutoff, isotope filter, peak annotations) and the data of the
        shown overview is used again.
        """
        if (settings := self._get_overview_settings()) is None:
            return
        key: tuple = (mode, tuple(f_ids), f_id_lib)
        data: OverviewData | None = self.overview_data if reuse_data and key == self.shown_overview else None
        match_result = None
        if mode == "match":
            match_result = next(m for m in self.matches[f_ids[0]] if m['feature_id'] == f_id_lib)
        annotations = {f_id: self._get_feature_annotation(f_id) for f_id in f_ids}

        self.match_panel_request += 1
        self.shown_overview = key
        self.pending_overview = (self.match_panel_request, key, match_result, annotations, data, settings)
        if not reuse_data:
            self.tk_status_field.config(text=f"Plotting {mode} of feature(s) {', '.join(map(str, f_ids))} ...")
        if not self.overview_running:
            self._start_overview_worker()

    def redraw_overview(self):
        """Redraw the shown overview with the current settings."""
        if self.shown_overview is not None and self.overview_data is not None:
            self.request_overview(*self.shown_overview, reuse_data=True)

    def _start_overview_worker(self):
        request, key, match_result, annotations, data, settings = self.pending_overview
        self.pending_overview = None
        self.overview_running = True
        threading.Thread(
            target=self._overview_worker,
            args=(self.meas_db, self.lib_db, self.ms2_spectra, request, key, match_result, annotations, data,
                  settings),
            daemon=True
        ).start()

    def _overview_worker(self, meas_db: FeatureManagerDB, lib_db: Library, ms2_spectra: dict | None, request: int,
                         key: tuple, match_result: dict | None, annotations: dict, data: OverviewData | None,
                         settings: OverviewSettings):
        """Worker thread: gather the data (unless given) and build the figure (not attached to tk yet)."""
        fig, spectra, error_text = None, None, None
        try:
            if data is None:
                mode, f_ids, f_id_lib = key
                data = OverviewData(
                    mode=mode,
                    features=[gather_feature(meas_db, f_id, ms2_spectra, annotations[f_id]) for f_id in f_ids],
                    library=gather_library(lib_db, f_id_lib) if mode == "match" else None,
                    match_result=match_result,
                )
            fig = Figure(figsize=(12, 10))
            spectra = build_overview_figure(fig, data, settings)
        except Exception as error:
            error_text = f"{type(error).__name__}: {error}"

        try:
            self.tk_root.after(0, lambda: self._finish_overview(request, key, data, fig, spectra, error_text))
        except (tk.TclError, RuntimeError):  # window was closed in the meantime
            pass

    def _finish_overview(self, request: int, key: tuple, data: OverviewData | None, fig: Figure | None,
                         spectra: list[PlottedSpectrum] | None, error_text: str | None):
        """Main thread: show the figure built by the worker, unless newer content was requested in the meantime."""
        self.overview_running = False
        if self.pending_overview is not None:
            self._start_overview_worker()
            return
        if request != self.match_panel_request:
            return

        if error_text is not None:
            self.shown_overview = None  # requesting it again retries
            self.tk_status_field.config(text="Could not plot the overview.")
            messagebox.showerror("Plot Error", f"Could not plot the overview:\n{error_text}")
            return

        self.overview_data, self.overview_spectra = data, spectra
        mode, f_ids, f_id_lib = key
        title = {
            "match": f"Match: feature {f_ids[0]} -> library feature {f_id_lib}",
            "feature": f"Feature {f_ids[0]}",
            "comparison": f"Comparison: feature {f_ids[0]} vs {f_ids[-1]}",
        }[mode]
        self._show_in_match_panel(fig, title)
        status: str = f"Plotted {title[0].lower()}{title[1:]}."
        if self.tk_hide_isotope_peaks.get():
            status += f" Hid {self._count_hidden_isotope_peaks(data, spectra)} MS2 peak(s) within ±4 Da of the precursor."
        self.tk_status_field.config(text=status)

    @staticmethod
    def _count_hidden_isotope_peaks(data: OverviewData, spectra: list[PlottedSpectrum]) -> int:
        """Number of MS2 peaks removed by the isotope filter (all shown spectra)."""
        n_peaks: dict[SpectrumKey, int] = {("meas", f.f_id): len(f.ms2.mzs) for f in data.features if f.ms2 is not None}
        if data.library is not None and data.library.ms2 is not None:
            n_peaks[("lib", data.library.f_id)] = len(data.library.ms2.mzs)
        return sum(n_peaks[spectrum.key] - len(spectrum.mzs) for spectrum in spectra if spectrum.key in n_peaks)

    def on_click_overview(self, event):
        """Click on an MS2 peak: edit its annotation (empty annotation shows the m/z again)."""
        if self.match_toolbar.mode or event.button != 1 or event.inaxes is None or self.overview_spectra is None:
            return
        # upper or lower (mirror) spectrum, depending on where the click is
        candidates = [spectrum for spectrum in self.overview_spectra
                      if spectrum.ax is event.inaxes and (spectrum.heights[0] > 0) == (event.ydata >= 0)]
        if not candidates:
            return
        spectrum = candidates[0]
        # closest peak in screen coordinates that reaches up to the click (small neighbors of a large peak are often
        # only a pixel away), at most a few pixels away
        x_pixels = event.inaxes.transData.transform(np.column_stack([spectrum.mzs, np.zeros(len(spectrum.mzs))]))[:, 0]
        distances = np.abs(x_pixels - event.x)
        distances[np.abs(spectrum.heights) + 0.1 < abs(event.ydata)] = np.inf
        idx: int = int(np.argmin(distances))
        if distances[idx] > 6:
            return

        mz: float = float(spectrum.mzs[idx])
        source, f_id = spectrum.key
        if source == "meas":
            # measured spectra: the annotation is stored in the measured database (and in the loaded spectra)
            peaks = self._get_loaded_peaks(f_id, mz)
            current: str = next((peak.annotation for peak in peaks if peak.annotation), "")
            where = "stored in the measured database"
        else:
            # library spectra are not changed, the annotation is only kept in this session
            current = self.peak_annotations.get(spectrum.key, {}).get(peak_key(mz), "")
            where = "only for this session"
        text = simpledialog.askstring(
            "Peak annotation", f"Annotation of the peak at m/z {mz:.4f} "
                               f"({'library feature' if source == 'lib' else 'feature'} {f_id}, {where}),\n"
                               f"empty to show the m/z again:",
            initialvalue=current, parent=self.tk_root
        )
        if text is None:  # cancelled
            return
        annotation: str | None = text.strip() or None
        if source == "meas":
            try:
                self.meas_db.set_peak_annotation(f_id, 2, mz, annotation)
            except Exception as error:
                messagebox.showerror("Peak annotation", f"Could not store the annotation:\n{error}")
                return
            for peak in peaks:
                peak.annotation = annotation
            self._update_fragment_choices()
        elif annotation:
            self.peak_annotations.setdefault(spectrum.key, {})[peak_key(mz)] = annotation
        else:
            self.peak_annotations.get(spectrum.key, {}).pop(peak_key(mz), None)
        self.redraw_overview()

    def _get_loaded_peaks(self, f_id: int, mz: float) -> list:
        """Peaks at mz of the MS2 spectra of a measured feature that are loaded (overview and cache of all spectra)."""
        spectra = [feature.ms2 for feature in (self.overview_data.features if self.overview_data else [])
                   if feature.f_id == f_id]
        if self.ms2_spectra is not None:
            spectra.append(self.ms2_spectra.get(f_id))
        peaks: dict[int, object] = {}  # by id, the overview usually uses the spectrum of the cache
        for spectrum in spectra:
            for peak in (spectrum.peaks if spectrum is not None else []):
                if peak_key(peak.mz) == peak_key(mz):
                    peaks[id(peak)] = peak
        return list(peaks.values())

    def on_change_annotation_cutoff(self, value: str):
        """Slider moved: update the label right away, redraw once the slider rests."""
        self.tk_annotation_cutoff_label.config(text=f"{float(value):.2f}")
        if self.annotation_cutoff_job is not None:
            self.tk_root.after_cancel(self.annotation_cutoff_job)
        self.annotation_cutoff_job = self.tk_root.after(300, self._apply_annotation_cutoff)

    def _apply_annotation_cutoff(self):
        self.annotation_cutoff_job = None
        self.redraw_overview()

    def show_feature(self, f_id: int):
        """Overview of a single measured feature (with or without matches)."""
        self.mark_selected_feature(f_id)
        self.request_overview("feature", [f_id])

    def on_submit_feature_ids(self):
        """Feature ID field: one ID shows that feature, two IDs (separated by , ; or space) compare them."""
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        try:
            f_ids = [int(part) for part in re.split(r"[,;\s]+", self.tk_feature_ids_entry.get().strip()) if part]
        except ValueError:
            messagebox.showerror("Invalid input", "Feature IDs must be whole numbers.")
            return
        if not 1 <= len(f_ids) <= 2 or len(set(f_ids)) != len(f_ids):
            messagebox.showerror("Invalid input", "Enter one feature ID or two different ones.")
            return
        if missing := [f_id for f_id in f_ids if f_id not in self.meas_db.mzs]:
            messagebox.showerror("Invalid input", f"There is no measured feature {', '.join(map(str, missing))}.")
            return
        if len(f_ids) == 1:
            self.clear_comparison()
            self.show_feature(f_ids[0])
        else:
            self.compared_f_ids = f_ids
            self.mark_compared_features()
            self.request_overview("comparison", f_ids)

    # ------------------------------------------------------------------------------------------------------------------
    # Scatter plot
    # ------------------------------------------------------------------------------------------------------------------
    def _get_rt_mz(self, f_id: int) -> tuple[float, float] | None:
        """RT in minutes (as in FeatureManagerDB.plot_scatter) and m/z of a measured feature."""
        rt = self.meas_db.retention_times_in_seconds.get(f_id)
        mz = self.meas_db.mzs.get(f_id)
        if rt is None or mz is None:
            return None
        return rt / 60, mz

    def _destroy_scatter_plot(self):
        if self.scatter_canvas is not None:
            self.scatter_canvas.get_tk_widget().destroy()
            self.scatter_canvas = None
        if self.scatter_toolbar is not None:
            self.scatter_toolbar.destroy()
            self.scatter_toolbar = None
        self.scatter_ax = None
        self.filtered_features_marker = None
        self.selected_feature_marker = None
        self.compared_features_markers = []
        self.mass_difference_lines = []
        self.annotated_features_markers = []

    def draw_scatter_plot(self, keep_view: bool = True):
        """
        (Re)draw the scatter plot of the measurement including annotations, highlights and markers.
        With keep_view, the zoomed / panned region of the previous plot is kept (not for a new database).
        """
        view = None
        if keep_view and self.scatter_ax is not None:
            view = self.scatter_ax.get_xlim(), self.scatter_ax.get_ylim()
        self._destroy_scatter_plot()

        fig = Figure(figsize=(7, 6))
        self.scatter_ax = fig.add_subplot(111)
        try:
            self.meas_db.plot_scatter(
                add_annotations=self.tk_show_annotations.get(),
                highlight_features=self.highlighted_f_ids,
                ax=self.scatter_ax
            )
        except Exception as error:
            self.scatter_ax = None
            messagebox.showerror("Error", f"Could not render the scatter plot:\n{error}")
            return
        self.scatter_ax.set_xlabel("Retention time (min)")
        self.scatter_ax.set_ylabel("m/z")
        fig.tight_layout()

        self.scatter_canvas = FigureCanvasTkAgg(fig, master=self.tk_scatter_plot_frame)
        self.scatter_toolbar = NavigationToolbar2Tk(self.scatter_canvas, self.tk_scatter_plot_frame)
        self.scatter_toolbar.update()
        self.scatter_toolbar.pack()
        self._pack_canvas(self.scatter_canvas)
        self.scatter_canvas.mpl_connect('button_press_event', self.on_click_scatter)
        self.scatter_canvas.mpl_connect('key_press_event', self.on_key_scatter)
        self.scatter_canvas.mpl_connect('figure_enter_event', self.on_enter_scatter)
        self.scatter_canvas.mpl_connect('figure_leave_event', self.on_leave_scatter)

        self.mark_annotated_features()
        self.mark_filtered_features()
        if (match := self._get_selected_match()) is not None:
            self.mark_selected_feature(int(match.feature_id_meas))
        self.mark_compared_features()
        self.mark_mass_differences()
        if view is not None:
            # the full view first, so that the home button of the toolbar still shows all features
            self.scatter_toolbar.push_current()
            self.scatter_ax.set_xlim(*view[0])
            self.scatter_ax.set_ylim(*view[1])
            self.scatter_toolbar.push_current()
        self.scatter_canvas.draw()

    def redraw_scatter_plot(self):
        """Redraw the scatter plot, if there is one (i.e. a measured database is loaded)."""
        if self.scatter_ax is not None:
            self.draw_scatter_plot()

    def on_click_highlight(self):
        """
        Highlight features whose MS2 spectrum contains the entered fragments and neutral losses.

        Peaks are matched with the MS2 fragment tolerance, peaks below the minimum peak intensity are ignored.
        With "all", a feature needs every entered fragment and loss, with "any" one of them is enough.
        Empty fields remove the highlight.
        """
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        try:
            fragments = parse_masses(self.tk_fragments_entry.get())
            losses = parse_masses(self.tk_losses_entry.get())
            mz_tol_da: float = self.tk_max_dmz_ms2_mda.get() * 1e-3
            min_intensity: float = self.tk_min_peak_intensity.get()
        except ValueError as error:
            messagebox.showerror("Invalid input", str(error))
            return
        except tk.TclError:
            messagebox.showerror("Invalid input", "Check the MS2 fragment tolerance and the minimum peak intensity.")
            return
        if not (fragments or losses):
            self.on_click_remove_highlight()
            return

        try:
            f_id_sets: list[set[int]] = (
                    [set(self.meas_db.filter_by_fragment_mz(mass, mz_tol_da, min_intensity)) for _, mass in fragments]
                    + [set(self.meas_db.filter_by_neutral_loss(mass, mz_tol_da, min_intensity)) for _, mass in losses]
            )
        except Exception as error:
            messagebox.showerror("Error", f"Could not search for fragments and losses:\n{error}")
            return
        mode: str = self.tk_highlight_mode.get()
        combined = set.intersection(*f_id_sets) if mode == "all" else set.union(*f_id_sets)
        self.highlighted_f_ids = sorted(combined)
        self.draw_scatter_plot()

        def describe(entries: list[tuple[str, float]]) -> str:
            return ", ".join(entry if entry == f"{mass:g}" else f"{entry} ({mass:.4f})"
                             for entry, mass in entries)

        searched = "; ".join(part for part in (
            f"fragments {describe(fragments)}" if fragments else "",
            f"losses {describe(losses)}" if losses else "",
        ) if part)
        self.tk_status_field.config(
            text=f"Highlighted {len(self.highlighted_f_ids)} features with {mode} of {searched} "
                 f"(peaks above {min_intensity:g})."
        )

    def _update_fragment_choices(self):
        """Fill the fragment list with the annotated MS2 peaks of the measured database, most frequent first."""
        self.fragment_choices = {}
        try:
            annotations = self.meas_db.get_peak_annotations(level=2) if self.meas_db is not None else []
        except Exception:  # e.g. not an msIO database
            annotations = []
        for annotation, mz, n_features in annotations:
            try:  # formulas (e.g. C5H15NO4P+) are inserted as they are, other names as their m/z
                parse_masses(annotation)
                value = annotation
            except ValueError:
                value = f"{mz:.4f}"
            label = f"{annotation}  (m/z {mz:.4f}, {n_features} feature{'s' if n_features > 1 else ''})"
            self.fragment_choices[label] = value
        self.tk_fragments_entry.config(values=list(self.fragment_choices))

    def _remember_fragments_text(self):
        # the combobox replaces the text with the selected entry, keep what was typed before
        self.fragments_text_before_selection = self.tk_fragments_entry.get()

    def on_select_fragment(self, event=None):
        """Append the selected annotated fragment to the fragments typed before."""
        value: str = self.fragment_choices.get(self.tk_fragments_entry.get(), self.tk_fragments_entry.get())
        before: str = self.fragments_text_before_selection.strip().rstrip(";").strip()
        self.tk_fragments_entry.set(f"{before}; {value}" if before else value)
        self.tk_fragments_entry.icursor("end")

    def on_click_remove_highlight(self):
        self.highlighted_f_ids = None
        self.redraw_scatter_plot()

    def on_click_mark_mass_differences(self):
        """
        Connect pairs of features whose m/z differ by the entered mass differences or adduct exchanges.

        m/z values are matched with the precursor tolerance (stricter of ppm and mDa), the retention times of
        both features have to be within the RT window. Empty field removes the marks.
        """
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        try:
            differences = parse_mass_differences(self.tk_mass_differences_entry.get())
            max_dmz_ppm: float = self.tk_max_dmz_ppm.get()
            max_dmz_da: float = self.tk_max_dmz_mda.get() * 1e-3
            rt_window: float = self.tk_rt_window.get()
        except ValueError as error:
            messagebox.showerror("Invalid input", str(error))
            return
        except tk.TclError:
            messagebox.showerror("Invalid input", "Check the precursor tolerances and the RT window.")
            return
        if not differences:
            self.on_click_remove_mass_differences()
            return

        points = {f_id: p for f_id in self.meas_db.feature_ids if (p := self._get_rt_mz(f_id)) is not None}
        f_ids = np.array(list(points.keys()))
        rts, mzs = np.array(list(points.values())).T
        self.mass_difference_pairs = [
            (entry, find_feature_pairs(f_ids, rts, mzs, expected_partner_mz, max_dmz_ppm, max_dmz_da, rt_window))
            for entry, expected_partner_mz in differences
        ]
        self.mark_mass_differences()

        counts = ", ".join(f"{entry}: {len(pairs)}" for entry, pairs in self.mass_difference_pairs)
        self.tk_status_field.config(text=f"Pairs of features within {rt_window:g} min: {counts}.")

    def on_click_remove_mass_differences(self):
        self.mass_difference_pairs = []
        self.mark_mass_differences()

    def find_marked_mass_differences(self, f_id_a: int, f_id_b: int) -> list[str]:
        """Marked mass differences that connect two features (in either direction)."""
        return [entry for entry, pairs in self.mass_difference_pairs
                if (f_id_a, f_id_b) in pairs or (f_id_b, f_id_a) in pairs]

    def mark_mass_differences(self):
        """Draw a line between the features of each pair with a marked mass difference."""
        if self.scatter_ax is None:
            return
        for artist in self.mass_difference_lines:
            artist.remove()
        self.mass_difference_lines = []

        for i, (entry, pairs) in enumerate(self.mass_difference_pairs):
            color: str = MASS_DIFFERENCE_COLORS[i % len(MASS_DIFFERENCE_COLORS)]
            segments = [(self._get_rt_mz(f_id_a), self._get_rt_mz(f_id_b)) for f_id_a, f_id_b in pairs]
            lines = LineCollection(segments, colors=color, linewidths=2.5, zorder=12, label=f"{entry} ({len(pairs)})")
            self.scatter_ax.add_collection(lines)
            self.mass_difference_lines.append(lines)
            # mass differences are small compared to the m/z range, the end caps keep short lines visible
            if segments:
                ends = np.array([point for segment in segments for point in segment])
                self.mass_difference_lines.append(self.scatter_ax.scatter(
                    ends[:, 0], ends[:, 1], marker="_", s=120, linewidths=2.5, color=color, zorder=12
                ))
        self._update_scatter_legend()

    def mark_annotated_features(self):
        """
        Color features with manual annotations (in the database or to add) green, highlighted ones in yellow
        (red + green), independent of whether the annotation labels are shown.
        """
        if self.scatter_ax is None:
            return
        for artist in self.annotated_features_markers:
            artist.remove()
        self.annotated_features_markers = []

        annotated: set[int] = {key[1] for collection in (self.saved_annotations, self.pending_annotations)
                               for key in collection if key not in self.discarded_annotations}
        highlighted: set[int] = set(self.highlighted_f_ids or [])
        for f_ids, color, label in ((annotated - highlighted, "green", "Annotated"),
                                    (annotated & highlighted, "gold", "Annotated + highlighted")):
            points = [p for f_id in f_ids if (p := self._get_rt_mz(f_id)) is not None]
            if points:
                rts, mzs = zip(*points)
                # same size as the highlighted features of FeatureManagerDB.plot_scatter
                self.annotated_features_markers.append(
                    self.scatter_ax.scatter(rts, mzs, c=color, s=20, zorder=3, label=label)
                )
        self._update_scatter_legend()

    def mark_filtered_features(self):
        """Circle features of the filtered table in black (only if a filter removes rows and circles are enabled)."""
        if self.scatter_ax is None:
            return
        if self.filtered_features_marker is not None:
            self.filtered_features_marker.remove()
        self.filtered_features_marker = None

        if self.tk_circle_matches.get() and 0 < len(self.filtered_matches_table) < len(self.matches_table):
            points = [self._get_rt_mz(f_id) for f_id in self.filtered_matches_table.feature_id_meas.unique()]
            points = [p for p in points if p is not None]
            if points:
                rts, mzs = zip(*points)
                self.filtered_features_marker = self.scatter_ax.scatter(
                    rts, mzs, s=120, facecolors='none', edgecolors='black', linewidths=2, zorder=10,
                    label='Library match'
                )
        self._update_scatter_legend()

    def mark_selected_feature(self, f_id_meas: int):
        """Circle the given feature in red."""
        if self.scatter_ax is None or (point := self._get_rt_mz(f_id_meas)) is None:
            return
        if self.selected_feature_marker is not None:
            self.selected_feature_marker.remove()

        rt, mz = point
        self.selected_feature_marker = self.scatter_ax.scatter(
            [rt], [mz], s=300, facecolors="none", edgecolors="red", linewidths=2, zorder=10, label='Selected feature'
        )
        self._update_scatter_legend()

    def mark_compared_features(self):
        """Put squares around the features selected by shift-click."""
        if self.scatter_ax is None:
            return
        for artist in self.compared_features_markers:
            artist.remove()
        self.compared_features_markers = []

        for f_id, color in zip(self.compared_f_ids, COMPARISON_COLORS):
            if (point := self._get_rt_mz(f_id)) is None:
                continue
            self.compared_features_markers.append(self.scatter_ax.scatter(
                [point[0]], [point[1]], s=200, marker="s", facecolors="none", edgecolors=color, linewidths=2,
                zorder=11, label=f"Compared feature {f_id}"
            ))
        self._update_scatter_legend()

    def _update_scatter_legend(self):
        # only the markers have labels, so the legend lists exactly the markers that are shown
        if self.scatter_ax.get_legend_handles_labels()[0]:
            self.scatter_ax.legend(loc='upper right')
        elif (legend := self.scatter_ax.get_legend()) is not None:
            legend.remove()
        self.scatter_canvas.draw_idle()

    def _find_closest_feature(self, x: float, y: float) -> int | None:
        """
        Feature closest to a point of the scatter plot (highlighted or not).

        RT and m/z have different scales, so the distance is measured relative to the visible axis range.
        """
        points = {f_id: p for f_id in self.meas_db.feature_ids if (p := self._get_rt_mz(f_id)) is not None}
        if not points:
            return None
        f_ids = np.array(list(points.keys()))
        rts, mzs = np.array(list(points.values())).T

        x_min, x_max = self.scatter_ax.get_xlim()
        y_min, y_max = self.scatter_ax.get_ylim()
        distances = ((rts - x) / (x_max - x_min)) ** 2 + ((mzs - y) / (y_max - y_min)) ** 2
        return int(f_ids[np.argmin(distances)])

    def on_click_scatter(self, event):
        """
        Click: select the best match of the closest feature in the results table and show it, features without
        matches (in the filtered table) are shown alone.
        Shift-click: add the closest feature to the comparison (the oldest of two is replaced).
        """
        # ignore clicks outside the axes or while zooming / panning with the toolbar
        if event.inaxes is not self.scatter_ax or self.scatter_toolbar.mode:
            return
        if (closest_f_id := self._find_closest_feature(event.xdata, event.ydata)) is None:
            return
        self.last_clicked_f_id = closest_f_id

        # modifiers come from the state of this click; event.key is the last key seen by the canvas and can get
        # stuck at "shift" when the key is released while another widget has the focus
        if "shift" in event.modifiers:
            self.add_feature_to_comparison(closest_f_id)
            return

        self.clear_comparison()
        # the overview is requested here, the selection event isn't generated if the row was already selected
        if self.select_best_match_in_table(closest_f_id, show_match=False):
            match = self._get_selected_match()
            self.mark_selected_feature(closest_f_id)
            self.request_overview("match", [closest_f_id], int(match.feature_id_lib))
        else:
            self.show_feature(closest_f_id)

    def select_best_match_in_table(self, f_id_meas: int, show_match: bool = True) -> bool:
        """
        Select the match with the highest MS2 score of a measured feature in the results table (selecting a row
        shows the match overview, unless show_match is False). Returns False if the filtered table has no match.
        """
        rows = self.filtered_matches_table.index[self.filtered_matches_table.feature_id_meas == f_id_meas]
        if rows.empty:
            return False
        # matches without MS2 score come last
        best_row = self.filtered_matches_table.loc[rows, 'ms2_score'].fillna(-np.inf).idxmax()
        # rows of the table are in the order of filtered_matches_table (its index is reset when filtering)
        item: str = self.tk_matches_table.get_children()[best_row]
        if not show_match:
            self._ignore_pending_table_select_events()
        self.tk_matches_table.selection_set(item)
        self.tk_matches_table.see(item)
        return True

    def on_enter_scatter(self, event):
        """
        Keys only reach the canvas if it has the focus: take it while the mouse is over the plot and give it back to
        the previous widget (e.g. a filter field) when the mouse leaves.
        """
        try:
            focused = self.tk_root.focus_get()
        except KeyError:  # tkinter can't resolve some internal widgets, e.g. the open list of a combobox
            focused = None
        canvas_widget = self.scatter_canvas.get_tk_widget()
        if focused is not canvas_widget:
            self.focus_before_scatter = focused
            canvas_widget.focus_set()

    def on_leave_scatter(self, event):
        canvas_widget = self.scatter_canvas.get_tk_widget()
        if self.focus_before_scatter is not None and self.tk_root.focus_get() is canvas_widget:
            try:
                self.focus_before_scatter.focus_set()
            except tk.TclError:  # widget was destroyed in the meantime
                pass
        self.focus_before_scatter = None

    def on_key_scatter(self, event):
        """+ / - zoom in / out by scatter_zoom_factor, centered on the mouse position."""
        # key names depend on the keyboard and on whether the key produces a character (e.g. numpad)
        if event.inaxes is not self.scatter_ax or event.key not in ZOOM_IN_KEYS | ZOOM_OUT_KEYS:
            return
        scale: float = 1 / self.scatter_zoom_factor if event.key in ZOOM_IN_KEYS else self.scatter_zoom_factor
        # the point under the mouse stays where it is
        x_min, x_max = self.scatter_ax.get_xlim()
        y_min, y_max = self.scatter_ax.get_ylim()
        self.scatter_ax.set_xlim(event.xdata - (event.xdata - x_min) * scale, event.xdata + (x_max - event.xdata) * scale)
        self.scatter_ax.set_ylim(event.ydata - (event.ydata - y_min) * scale, event.ydata + (y_max - event.ydata) * scale)
        # makes the zoom steps available to the back / forward buttons of the toolbar
        self.scatter_toolbar.push_current()
        self.scatter_canvas.draw_idle()

    def clear_comparison(self):
        """The overview shows something else than the shift-clicked features: remove them and their squares."""
        if self.compared_f_ids:
            self.compared_f_ids = []
            self.mark_compared_features()

    def add_feature_to_comparison(self, f_id: int):
        """Add a feature to the comparison: the first one is shown alone, two are compared."""
        if f_id in self.compared_f_ids:
            self.tk_status_field.config(text=f"Feature {f_id} is already compared.")
            return
        self.compared_f_ids = (self.compared_f_ids + [f_id])[-2:]
        self.mark_compared_features()
        if len(self.compared_f_ids) < 2:
            self.request_overview("feature", [f_id])
        else:
            self.request_overview("comparison", list(self.compared_f_ids))


class AnnotationManagerWindow(tk.Toplevel):
    """
    Window to review manual annotations before writing them to the measured database.

    Lists new annotations (library matches, manual annotations, in-source fragments) and the ones
    already in the database. Each can be kept or discarded; "Write to database" adds the kept new
    ones and removes the discarded ones from the database. The state lives in the AnnotationGUI,
    so closing this window doesn't lose any changes.
    """

    def __init__(self, gui: AnnotationGUI):
        super().__init__(gui.tk_root)
        self.gui: AnnotationGUI = gui
        self.title("Manual annotations")
        self.geometry("1000x500")
        self.protocol("WM_DELETE_WINDOW", self.on_close)
        # rows are identified by their position, refresh() rebuilds this mapping
        self.iid_to_key: dict[str, AnnotationKey] = {}

        ttk.Label(
            self,
            text="Double-click a row or press Space to keep / discard it. \"Write to database\" adds the kept new "
                 "annotations and removes the discarded ones from the database.",
            wraplength=960, padding=(10, 10, 10, 5)
        ).pack(fill="x")

        tk_table_frame = ttk.Frame(self, padding=(10, 0))
        tk_table_frame.pack(fill="both", expand=True)
        columns: dict[str, tuple[str, int, str]] = {
            "keep": ("Keep", 45, "center"),
            "state": ("State", 85, "center"),
            "kind": ("Kind", 110, "w"),
            "feature_id_meas": ("Meas. ID", 70, "center"),
            "reference": ("Reference", 90, "center"),
            "name": ("Name", 250, "w"),
            "formula": ("Formula", 110, "w"),
            "adduct": ("Adduct", 85, "w"),
            "ms2_score": ("MS2 score", 75, "e"),
        }
        self.tk_table = ttk.Treeview(tk_table_frame, columns=list(columns), show="headings", selectmode="extended")
        for column, (heading, width, anchor) in columns.items():
            self.tk_table.heading(column, text=heading)
            self.tk_table.column(column, width=width, anchor=anchor, stretch=column == "name")
        self.tk_table.tag_configure("to add", background="#CCFFCC")
        self.tk_table.tag_configure("to remove", background="#FFCCCC")
        self.tk_table.tag_configure("discarded", foreground="grey")
        self.tk_table.bind("<Double-1>", self.on_double_click)
        self.tk_table.bind("<space>", lambda event: self.gui.toggle_annotations(self._get_selected_keys()))

        y_scroll = ttk.Scrollbar(tk_table_frame, orient="vertical", command=self.tk_table.yview)
        self.tk_table.configure(yscrollcommand=y_scroll.set)
        self.tk_table.pack(side="left", fill="both", expand=True)
        y_scroll.pack(side="left", fill="y")

        tk_buttons = ttk.Frame(self, padding=10)
        tk_buttons.pack(fill="x")
        ttk.Button(tk_buttons, text="Keep selected",
                   command=lambda: self.gui.toggle_annotations(self._get_selected_keys(), keep=True)).pack(side="left")
        ttk.Button(tk_buttons, text="Discard selected",
                   command=lambda: self.gui.toggle_annotations(self._get_selected_keys(), keep=False)).pack(
            side="left", padx=5)
        ttk.Button(tk_buttons, text="Keep all",
                   command=lambda: self.gui.toggle_annotations(list(self.iid_to_key.values()), keep=True)).pack(
            side="left")
        ttk.Button(tk_buttons, text="Discard all",
                   command=lambda: self.gui.toggle_annotations(list(self.iid_to_key.values()), keep=False)).pack(
            side="left", padx=5)
        ttk.Button(tk_buttons, text="Close", command=self.on_close).pack(side="right")
        ttk.Button(tk_buttons, text="Write to database",
                   command=lambda: self.gui.write_annotations_to_db(parent=self)).pack(side="right", padx=5)

        self.refresh()

    def _get_selected_keys(self) -> list[AnnotationKey]:
        return [self.iid_to_key[iid] for iid in self.tk_table.selection()]

    def refresh(self):
        """Rebuild the rows from the annotation state of the GUI, keeping selection and scroll position."""
        selected_keys = set(self._get_selected_keys())
        y_view = self.tk_table.yview()[0]
        self.tk_table.delete(*self.tk_table.get_children())
        self.iid_to_key = {}

        # new annotations first, then the ones in the database, each sorted by measured feature
        def sort_key(item) -> tuple:
            (kind, f_id_meas, ref), _ = item
            return f_id_meas, kind, ref

        rows = ([(key, info, False) for key, info in sorted(self.gui.pending_annotations.items(), key=sort_key)]
                + [(key, info, True) for key, info in sorted(self.gui.saved_annotations.items(), key=sort_key)])
        for idx, (key, info, is_saved) in enumerate(rows):
            kind, f_id_meas, _ = key
            discarded: bool = key in self.gui.discarded_annotations
            if is_saved:
                state = "to remove" if discarded else "in database"
            else:
                state = "discarded" if discarded else "to add"
            reference = info.get("reference")
            ms2_score = info.get("ms2_score")
            iid = str(idx)
            self.iid_to_key[iid] = key
            self.tk_table.insert("", "end", iid=iid, tags=(state,), values=(
                "✗" if discarded else "✓",
                state,
                ANNOTATION_KIND_LABELS.get(kind, kind),
                f_id_meas,
                "" if reference is None else f"{'lib.' if kind == 'library' else 'parent'} {reference}",
                info.get("name") or "",
                info.get("formula") or "",
                info.get("adduct") or "",
                "" if ms2_score is None or pd.isna(ms2_score) else f"{ms2_score:.3f}",
            ))

        self.tk_table.selection_set([iid for iid, key in self.iid_to_key.items() if key in selected_keys])
        self.tk_table.yview_moveto(y_view)

    def on_double_click(self, event):
        if iid := self.tk_table.identify_row(event.y):
            self.gui.toggle_annotations([self.iid_to_key[iid]])

    def on_close(self):
        self.gui.annotation_window = None
        self.destroy()


class ManualAnnotationDialog(tk.Toplevel):
    """
    Dialog to annotate a measured feature with a compound entered by hand.

    Name, formula, adduct and SMILES can be typed in or taken from the IPL builder of LipidCalculator.
    The m/z calculated from formula and adduct is compared to the m/z of the feature.
    """

    def __init__(self, gui: AnnotationGUI, f_id_meas: int | None):
        super().__init__(gui.tk_root)
        self.gui: AnnotationGUI = gui
        self.title("Add manual annotation")
        self.resizable(True, False)
        self.tk_f_id_meas = tk.StringVar(self, value="" if f_id_meas is None else str(f_id_meas))
        self.tk_name = tk.StringVar(self)
        self.tk_formula = tk.StringVar(self)
        self.tk_adduct = tk.StringVar(self, value="[M+H]+")
        self.tk_smiles = tk.StringVar(self)
        self.tk_mz_info = tk.StringVar(self)

        tk_form = ttk.Frame(self, padding=10)
        tk_form.pack(fill="x")
        rows = [("Measured feature ID", ttk.Entry(tk_form, textvariable=self.tk_f_id_meas, width=12)),
                ("Name", ttk.Entry(tk_form, textvariable=self.tk_name, width=50)),
                ("Formula (neutral)", ttk.Entry(tk_form, textvariable=self.tk_formula, width=30)),
                ("Adduct", ttk.Combobox(tk_form, textvariable=self.tk_adduct, values=ADDUCTS, width=15)),
                ("SMILES (optional)", ttk.Entry(tk_form, textvariable=self.tk_smiles, width=50))]
        for row, (label, widget) in enumerate(rows):
            ttk.Label(tk_form, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            widget.grid(row=row, column=1, sticky="w")
        ttk.Label(tk_form, textvariable=self.tk_mz_info, foreground="gray").grid(
            row=len(rows), column=0, columnspan=2, sticky="w", pady=(5, 0))
        for variable in (self.tk_f_id_meas, self.tk_formula, self.tk_adduct):
            variable.trace_add("write", lambda *_: self._update_mz_info())

        tk_buttons = ttk.Frame(self, padding=(10, 0, 10, 10))
        tk_buttons.pack(fill="x")
        ttk.Button(tk_buttons, text="Build IPL ...", command=self.on_click_build_ipl).pack(side="left")
        ttk.Button(tk_buttons, text="Cancel", command=self.destroy).pack(side="right")
        ttk.Button(tk_buttons, text="Add", command=self.on_click_add).pack(side="right", padx=5)
        self._update_mz_info()

    def _update_mz_info(self):
        """Show the m/z calculated from formula and adduct next to the measured m/z of the feature."""
        parts: list[str] = []
        mz_meas = None
        try:
            mz_meas = self.gui.meas_db.mzs.get(int(self.tk_f_id_meas.get()))
            if mz_meas is not None:
                parts.append(f"measured m/z {mz_meas:.4f}")
        except ValueError:
            pass
        try:
            mz_calc = Adduct(self.tk_adduct.get().strip()).mass_to_mz(CompoundDict(self.tk_formula.get().strip()).mass)
            parts.append(f"calculated m/z {mz_calc:.4f}")
            if mz_meas is not None:
                parts.append(f"Δ {(mz_meas - mz_calc) / mz_calc * 1e6:.1f} ppm")
        except Exception:  # formula or adduct incomplete
            pass
        self.tk_mz_info.set(", ".join(parts))

    def on_click_build_ipl(self):
        """Open the IPL builder of LipidCalculator, its "Use for annotation" button fills this dialog."""
        open_ipl_builder(self, "Use for annotation", self._use_ipl)

    def _use_ipl(self, frame, payload: dict, tk_window: tk.Toplevel):
        """Take name, formula, SMILES and the selected adduct (if any) from the IPL builder."""
        self.tk_name.set(frame.name_var.get())
        self.tk_formula.set(payload["formula"])
        if payload.get("adduct"):
            self.tk_adduct.set(payload["adduct"])
        self.tk_smiles.set(frame.smiles_var.get())
        tk_window.destroy()
        self.lift()

    def on_click_add(self):
        try:
            f_id_meas = int(self.tk_f_id_meas.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Measured feature ID must be a whole number.", parent=self)
            return
        if f_id_meas not in self.gui.meas_db.mzs:
            messagebox.showerror("Invalid input", f"There is no measured feature {f_id_meas}.", parent=self)
            return
        name, formula = self.tk_name.get().strip(), self.tk_formula.get().strip() or None
        adduct, smiles = self.tk_adduct.get().strip() or None, self.tk_smiles.get().strip() or None
        if not name:
            messagebox.showerror("Invalid input", "Please enter a name.", parent=self)
            return
        try:
            if formula:
                CompoundDict(formula)
            if adduct:
                Adduct(adduct)
        except Exception as error:
            messagebox.showerror("Invalid input", f"Invalid formula or adduct:\n{error}", parent=self)
            return

        added: bool = self.gui.add_pending_annotation(("manual", f_id_meas, f"{name}|{adduct}"), dict(
            kind="manual", name=name, formula=formula, adduct=adduct, smiles=smiles, reference=None
        ))
        self.gui.tk_status_field.config(
            text=f"Added manual annotation {name} {adduct or ''} for feature {f_id_meas}." if added
            else f"Feature {f_id_meas} is already annotated as {name} {adduct or ''}."
        )
        self.destroy()


class MzFromCompoundDialog(tk.Toplevel):
    """
    Dialog to set the m/z filter of the results table from a compound: searched by name in the library, built with
    the IPL builder of LipidCalculator or entered as formula and adduct. The filter is m/z ± precursor tolerance.
    """

    MAX_SEARCH_RESULTS: int = 500

    def __init__(self, gui: AnnotationGUI):
        super().__init__(gui.tk_root)
        self.gui: AnnotationGUI = gui
        self.title("m/z filter from compound")
        self.geometry("800x560")
        self.tk_name = tk.StringVar(self)
        self.tk_exact = tk.BooleanVar(self, value=False)
        self.tk_search_info = tk.StringVar(self)
        self.tk_formula = tk.StringVar(self)
        self.tk_adduct = tk.StringVar(self, value="[M+H]+")
        self.tk_mz = tk.StringVar(self)
        self.tk_range_info = tk.StringVar(self)
        # library search results: {row id in the table: (formula, adduct, m/z of the library)}
        self.search_results: dict[str, tuple[str | None, str | None, float | None]] = {}

        tk_search = ttk.LabelFrame(self, text="Search the library by name", padding=8)
        tk_search.pack(fill="both", expand=True, padx=10, pady=(10, 5))
        tk_search_bar = ttk.Frame(tk_search)
        tk_search_bar.pack(fill="x")
        tk_name_entry = ttk.Entry(tk_search_bar, textvariable=self.tk_name, width=40)
        tk_name_entry.pack(side="left")
        tk_name_entry.bind("<Return>", lambda event: self.on_click_search())
        ttk.Checkbutton(tk_search_bar, text="exact name", variable=self.tk_exact).pack(side="left", padx=8)
        self.tk_search_button = ttk.Button(tk_search_bar, text="Search", command=self.on_click_search)
        self.tk_search_button.pack(side="left")
        ttk.Label(tk_search_bar, textvariable=self.tk_search_info, foreground="gray").pack(side="left", padx=8)

        tk_table_frame = ttk.Frame(tk_search)
        tk_table_frame.pack(fill="both", expand=True, pady=(5, 0))
        columns: dict[str, tuple[str, int, str]] = {
            "name": ("Name", 330, "w"), "formula": ("Formula", 140, "w"), "adduct": ("Adduct", 90, "w"),
            "mz": ("m/z", 90, "e"), "f_id": ("Lib. ID", 70, "center"),
        }
        self.tk_table = ttk.Treeview(tk_table_frame, columns=list(columns), show="headings", selectmode="browse",
                                     height=10)
        for column, (heading, width, anchor) in columns.items():
            self.tk_table.heading(column, text=heading)
            self.tk_table.column(column, width=width, anchor=anchor, stretch=column == "name")
        self.tk_table.bind("<<TreeviewSelect>>", self.on_select_result)
        y_scroll = ttk.Scrollbar(tk_table_frame, orient="vertical", command=self.tk_table.yview)
        self.tk_table.configure(yscrollcommand=y_scroll.set)
        self.tk_table.pack(side="left", fill="both", expand=True)
        y_scroll.pack(side="left", fill="y")

        tk_compound = ttk.LabelFrame(self, text="Compound (from the search, the IPL builder or typed in)", padding=8)
        tk_compound.pack(fill="x", padx=10, pady=5)
        ttk.Label(tk_compound, text="Formula (neutral)").grid(row=0, column=0, sticky="w")
        ttk.Entry(tk_compound, textvariable=self.tk_formula, width=24).grid(row=0, column=1, sticky="w", padx=(4, 12))
        ttk.Label(tk_compound, text="Adduct").grid(row=0, column=2, sticky="w")
        ttk.Combobox(tk_compound, textvariable=self.tk_adduct, values=ADDUCTS, width=14).grid(
            row=0, column=3, sticky="w", padx=(4, 12))
        ttk.Button(tk_compound, text="Build IPL ...", command=self.on_click_build_ipl).grid(row=0, column=4)
        ttk.Label(tk_compound, text="m/z").grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(tk_compound, textvariable=self.tk_mz, width=14).grid(row=1, column=1, sticky="w", padx=(4, 12),
                                                                       pady=(6, 0))
        ttk.Label(tk_compound, textvariable=self.tk_range_info, foreground="gray").grid(
            row=1, column=2, columnspan=3, sticky="w", pady=(6, 0))
        for variable in (self.tk_formula, self.tk_adduct):
            variable.trace_add("write", lambda *_: self._calculate_mz())
        self.tk_mz.trace_add("write", lambda *_: self._update_range_info())

        tk_buttons = ttk.Frame(self, padding=(10, 0, 10, 10))
        tk_buttons.pack(fill="x")
        ttk.Button(tk_buttons, text="Close", command=self.destroy).pack(side="right")
        ttk.Button(tk_buttons, text="Apply m/z filter", command=self.on_click_apply).pack(side="right", padx=5)

        if self.gui.lib_db is None:
            self.tk_search_button.state(["disabled"])
            self.tk_search_info.set("load a library database to search it")
        tk_name_entry.focus_set()

    def _calculate_mz(self) -> bool:
        """m/z from formula and adduct, False if they can't be parsed (the m/z field is left as it is)."""
        try:
            mz = Adduct(self.tk_adduct.get().strip()).mass_to_mz(CompoundDict(self.tk_formula.get().strip()).mass)
        except Exception:
            return False
        self.tk_mz.set(f"{mz:.4f}")
        return True

    def _update_range_info(self):
        try:
            mz = float(self.tk_mz.get())
            tolerance = self.gui.get_precursor_tolerance_da(mz)
        except (ValueError, tk.TclError):
            self.tk_range_info.set("")
            return
        self.tk_range_info.set(f"filter: {mz - tolerance:.4f} - {mz + tolerance:.4f} (± {tolerance * 1e3:.2f} mDa, "
                               f"precursor tolerance)")

    def on_click_search(self):
        """Search the library names in a worker thread (millions of names in large libraries)."""
        if not (name := self.tk_name.get().strip()) or self.gui.lib_db is None:
            return
        self.tk_search_button.state(["disabled"])
        self.tk_search_info.set("searching ...")
        lib_db, exact = self.gui.lib_db, self.tk_exact.get()

        def worker():
            try:
                f_ids = lib_db.find_by_name(name, substring=not exact)
                rows = [(f_id, lib_db.names.get(f_id), lib_db.formula_metaboscape.get(f_id), lib_db.adduct.get(f_id),
                         lib_db.mzs.get(f_id)) for f_id in f_ids[:self.MAX_SEARCH_RESULTS]]
                result = (len(f_ids), rows, None)
            except Exception as error:
                result = (0, [], f"{type(error).__name__}: {error}")
            self.gui._call_in_main_thread(lambda: self._show_search_results(*result))

        threading.Thread(target=worker, daemon=True).start()

    def _show_search_results(self, n_found: int, rows: list[tuple], error_text: str | None):
        if not self.winfo_exists():  # closed meanwhile
            return
        self.tk_search_button.state(["!disabled"])
        self.tk_table.delete(*self.tk_table.get_children())
        self.search_results = {}
        if error_text is not None:
            self.tk_search_info.set(f"search failed: {error_text}")
            return
        for f_id, name, formula, adduct, mz in rows:
            iid = self.tk_table.insert("", "end", values=(name or "", formula or "", adduct or "",
                                                          "" if mz is None else f"{mz:.4f}", f_id))
            self.search_results[iid] = (formula, adduct, mz)
        shown = f", showing the first {len(rows)}" if n_found > len(rows) else ""
        self.tk_search_info.set(f"{n_found} compound(s) found{shown}")

    def on_select_result(self, event=None):
        """Take formula and adduct of the selected compound, the library m/z if they can't be parsed."""
        if not (selection := self.tk_table.selection()):
            return
        formula, adduct, mz = self.search_results[selection[0]]
        self.tk_formula.set(formula or "")
        self.tk_adduct.set(adduct or "")
        if not self._calculate_mz() and mz is not None:
            self.tk_mz.set(f"{mz:.4f}")

    def on_click_build_ipl(self):
        open_ipl_builder(self, "Use for m/z filter", self._use_ipl)

    def _use_ipl(self, frame, payload: dict, tk_window: tk.Toplevel):
        self.tk_formula.set(payload["formula"])
        if payload.get("adduct"):
            self.tk_adduct.set(payload["adduct"])
        tk_window.destroy()
        self.lift()

    def on_click_apply(self):
        try:
            mz = float(self.tk_mz.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Enter or select a compound with a valid m/z.", parent=self)
            return
        self.gui.set_mz_filter(mz)


class InSourceFragmentDialog(tk.Toplevel):
    """
    Dialog to annotate a measured feature as in-source fragment of another one.

    Shows the differences in m/z and RT of both features and the marked mass differences connecting them.
    """

    def __init__(self, gui: AnnotationGUI, f_id_fragment: int | None, f_id_parent: int | None):
        super().__init__(gui.tk_root)
        self.gui: AnnotationGUI = gui
        self.title("Mark in-source fragment")
        self.resizable(True, False)
        self.tk_fragment = tk.StringVar(self, value="" if f_id_fragment is None else str(f_id_fragment))
        self.tk_parent = tk.StringVar(self, value="" if f_id_parent is None else str(f_id_parent))
        self.tk_note = tk.StringVar(self)
        self.tk_info = tk.StringVar(self)

        tk_form = ttk.Frame(self, padding=10)
        tk_form.pack(fill="x")
        ttk.Label(tk_form, text="Fragment feature ID").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(tk_form, textvariable=self.tk_fragment, width=12).grid(row=0, column=1, sticky="w")
        ttk.Button(tk_form, text="Swap", command=self.on_click_swap).grid(row=0, column=2, rowspan=2, padx=8)
        ttk.Label(tk_form, text="Parent feature ID").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(tk_form, textvariable=self.tk_parent, width=12).grid(row=1, column=1, sticky="w")
        ttk.Label(tk_form, text="Note (optional)").grid(row=2, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(tk_form, textvariable=self.tk_note, width=40).grid(row=2, column=1, columnspan=2, sticky="w")
        ttk.Label(tk_form, textvariable=self.tk_info, foreground="gray", justify="left").grid(
            row=3, column=0, columnspan=3, sticky="w", pady=(5, 0))
        for variable in (self.tk_fragment, self.tk_parent):
            variable.trace_add("write", lambda *_: self._update_info())

        tk_buttons = ttk.Frame(self, padding=(10, 0, 10, 10))
        tk_buttons.pack(fill="x")
        ttk.Button(tk_buttons, text="Cancel", command=self.destroy).pack(side="right")
        ttk.Button(tk_buttons, text="Add", command=self.on_click_add).pack(side="right", padx=5)
        self._update_info()

    def _read_feature_ids(self) -> tuple[int, int] | None:
        try:
            f_ids = int(self.tk_fragment.get()), int(self.tk_parent.get())
        except ValueError:
            return None
        return f_ids if all(f_id in self.gui.meas_db.mzs for f_id in f_ids) else None

    def _update_info(self):
        """Show Δm/z, ΔRT and the marked mass differences connecting both features, suggest the note from them."""
        if (f_ids := self._read_feature_ids()) is None:
            self.tk_info.set("Enter two measured feature IDs (shift-click both in the scatter plot to prefill).")
            return
        (rt_fragment, mz_fragment), (rt_parent, mz_parent) = (self.gui._get_rt_mz(f_id) for f_id in f_ids)
        lines = [f"Δm/z (parent - fragment) {mz_parent - mz_fragment:.4f}, ΔRT {rt_parent - rt_fragment:+.2f} min"]
        if mz_parent < mz_fragment:
            lines.append("The parent has a lower m/z than the fragment.")
        if marked := self.gui.find_marked_mass_differences(*f_ids):
            lines.append(f"Marked mass differences: {', '.join(marked)}")
            if not self.tk_note.get():
                self.tk_note.set(f"loss of {', '.join(marked)}")
        self.tk_info.set("\n".join(lines))

    def on_click_swap(self):
        fragment, parent = self.tk_fragment.get(), self.tk_parent.get()
        self.tk_fragment.set(parent)
        self.tk_parent.set(fragment)

    def on_click_add(self):
        if (f_ids := self._read_feature_ids()) is None:
            messagebox.showerror("Invalid input", "Please enter two existing measured feature IDs.", parent=self)
            return
        f_id_fragment, f_id_parent = f_ids
        if f_id_fragment == f_id_parent:
            messagebox.showerror("Invalid input", "Fragment and parent must be different features.", parent=self)
            return
        note: str | None = self.tk_note.get().strip() or None
        added: bool = self.gui.add_pending_annotation(("in_source_fragment", f_id_fragment, str(f_id_parent)), dict(
            kind="in_source_fragment",
            name=f"in-source fragment of {f_id_parent}" + (f" ({note})" if note else ""),
            note=note,
            reference=f_id_parent,
        ))
        self.gui.tk_status_field.config(
            text=f"Marked feature {f_id_fragment} as in-source fragment of {f_id_parent}." if added
            else f"Feature {f_id_fragment} is already marked as in-source fragment of {f_id_parent}."
        )
        self.destroy()


if __name__ == "__main__":
    root: tk.Tk = tk.Tk()
    app: AnnotationGUI = AnnotationGUI(root)
    app.tk_root.mainloop()
