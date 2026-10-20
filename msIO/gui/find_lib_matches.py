"""
GUI for annotating measured features with matches from a spectral library.

Workflow:
    1. Select a measured database (FeatureManagerDB) and a library database (Library).
       The scatter plot (m/z vs RT) of the measurement is drawn right away and the MS2
       spectra of the measurement are loaded in the background.
    2. Set the matching parameters and click "Find Matches".
    3. Narrow down the results table with the filters (applied while typing) or by
       clicking a feature in the scatter plot.
    4. Select a row to plot the match overview, then "Annotate Feature" to store the
       library match as a manual annotation in the measured database.

Slow work (loading MS2 spectra, building the match overview figure) runs in worker
threads. Worker threads never touch tk widgets; they hand their results back to the
main thread with ``tk_root.after``.
"""

import numpy as np
import pandas as pd
import threading
import tkinter as tk
from tkinter import messagebox, filedialog, ttk

from LipidCalculator import CompoundDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from msIO.feature_managers.db import Library, FeatureManagerDB, lib_match_result_to_table

SQLITE_FILE_TYPES: list[tuple[str, str]] = [("SQLite Files", "*.sqlite"), ("All files", "*.*")]


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

        # database selection
        self.tk_meas_path: tk.StringVar = tk.StringVar()
        self.tk_lib_path: tk.StringVar = tk.StringVar()
        self.meas_db: FeatureManagerDB = None
        self.lib_db: Library = None
        # incremented for every loaded measured database, so results of workers for an older database are dropped
        self.meas_db_load_number: int = 0

        # MS2 spectra of the measurement, loaded in a worker thread
        self.ms2_spectra: dict = None
        self.ms2_loading: bool = False
        self.ms2_loading_error: str = None

        # matching parameters
        self.tk_max_dmz_ppm: tk.DoubleVar = tk.DoubleVar(value=10)
        self.tk_max_dmz_mda: tk.DoubleVar = tk.DoubleVar(value=3)
        self.tk_max_dmz_ms2_mda: tk.DoubleVar = tk.DoubleVar(value=10)
        self.tk_min_ms2_score: tk.DoubleVar = tk.DoubleVar(value=0)
        self.tk_require_ms2: tk.BooleanVar = tk.BooleanVar(value=False)

        # matching results: raw output of Library.find_matches and as table (all and filtered)
        self.matches: dict[int, list[dict]] = {}
        self.matches_table: pd.DataFrame = pd.DataFrame()
        self.filtered_matches_table: pd.DataFrame = pd.DataFrame()

        # filters are applied while typing, after a short pause
        self.live_filter_delay_ms: int = 300
        self.live_filter_job: str = None
        self.last_filter_values: tuple[str, ...] = None

        # match overview, the figure is built in a worker thread
        self.match_canvas: FigureCanvasTkAgg = None
        self.match_toolbar: NavigationToolbar2Tk = None
        self.match_plot_running: bool = False
        self.pending_match_plot: dict = None  # kwargs for Library.plot_match of the latest request

        # scatter plot (m/z vs RT) of the measurement
        self.scatter_canvas: FigureCanvasTkAgg = None
        self.scatter_toolbar: NavigationToolbar2Tk = None
        self.scatter_ax = None
        self.tk_show_annotations: tk.BooleanVar = tk.BooleanVar(value=True)
        self.neutral_loss_f_ids: list[int] = None  # features highlighted by the neutral loss filter
        self.filtered_features_marker = None  # black circles around features in the filtered table
        self.selected_feature_marker = None  # red circle around the feature selected in the table

        # background colors of table rows, alternating between measured features
        self.table_row_colors: list[str] = ["#E5FFCC", "#FFCCCC", "#CCFFCC", "#CCFFFF", "#CCCCFF", "#FFE5CC"]

        self._create_gui()

    # ------------------------------------------------------------------------------------------------------------------
    # GUI layout
    # ------------------------------------------------------------------------------------------------------------------
    def _create_gui(self):
        """
        Layout:
            top:    database selection | matching parameters
            middle: filters, results table, annotate button | match overview above scatter plot
            bottom: status bar
        """
        tk_top_frame = ttk.Frame(self.tk_root)
        tk_top_frame.pack(fill="x", padx=10, pady=(10, 5))
        self._create_database_selection(tk_top_frame)
        self._create_matching_parameters(tk_top_frame)

        # status bar is packed before the main area, so it keeps its space when the window shrinks
        self.tk_status_field = ttk.Label(
            self.tk_root, text="Please load both databases to begin.", relief="sunken", anchor="w", padding=(5, 2)
        )
        self.tk_status_field.pack(side="bottom", fill="x")

        tk_main_area = ttk.PanedWindow(self.tk_root, orient="horizontal")
        tk_main_area.pack(fill="both", expand=True, padx=10, pady=5)

        tk_results_pane = ttk.Frame(tk_main_area)
        tk_main_area.add(tk_results_pane, weight=1)
        self._create_filters(tk_results_pane)
        self._create_matches_table(tk_results_pane)
        ttk.Button(tk_results_pane, text="Annotate Feature", command=self.on_click_annotate).pack(
            anchor="e", pady=(5, 0)
        )

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
        tk_database_frame = ttk.LabelFrame(parent, text="Databases", padding=10)
        tk_database_frame.pack(side="left", fill="y", padx=(0, 10))

        ttk.Label(tk_database_frame, text="Measured database").grid(row=0, column=0, sticky="w")
        ttk.Entry(tk_database_frame, textvariable=self.tk_meas_path, width=45).grid(row=0, column=1, padx=5, pady=2)
        ttk.Button(tk_database_frame, text="Browse", command=self.on_browse_meas_db).grid(row=0, column=2)

        ttk.Label(tk_database_frame, text="Library database").grid(row=1, column=0, sticky="w")
        ttk.Entry(tk_database_frame, textvariable=self.tk_lib_path, width=45).grid(row=1, column=1, padx=5, pady=2)
        ttk.Button(tk_database_frame, text="Browse", command=self.on_browse_lib_db).grid(row=1, column=2)

    def _create_matching_parameters(self, parent):
        tk_parameter_frame = ttk.LabelFrame(parent, text="Matching parameters", padding=10)
        tk_parameter_frame.pack(side="left", fill="both", expand=True)

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

        self.tk_filter_min_f_id_meas, self.tk_filter_max_f_id_meas = add_range_filter(0, 0, "Meas. feature ID")
        self.tk_filter_min_rt, self.tk_filter_max_rt = add_range_filter(0, 1, "RT (min)")
        self.tk_filter_min_mz, self.tk_filter_max_mz = add_range_filter(1, 0, "m/z")
        self.tk_filter_subformula = add_filter(1, 1, "Formula contains")
        self.tk_filter_name = add_filter(2, 0, "Name contains (a; b)")
        self.tk_filter_min_ms2_score = add_filter(2, 1, "Min. MS2 score")
        self.tk_filter_min_n_hits_ms2 = add_filter(3, 0, "Min. #MS2 peak matches")

        tk_filter_buttons = ttk.Frame(tk_filter_frame)
        tk_filter_buttons.grid(row=3, column=4, columnspan=4, sticky="e")
        ttk.Button(tk_filter_buttons, text="Update filters", command=self.apply_filters).pack(side="left", padx=5)
        ttk.Button(tk_filter_buttons, text="Remove filters", command=self.on_click_remove_filters).pack(side="left")

        self.tk_filter_entries: list[ttk.Entry] = [
            self.tk_filter_min_f_id_meas, self.tk_filter_max_f_id_meas, self.tk_filter_name,
            self.tk_filter_min_rt, self.tk_filter_max_rt, self.tk_filter_subformula,
            self.tk_filter_min_ms2_score, self.tk_filter_min_n_hits_ms2, self.tk_filter_min_mz, self.tk_filter_max_mz
        ]
        for tk_entry in self.tk_filter_entries:
            tk_entry.bind("<KeyRelease>", self.on_edit_filter)

    def _create_matches_table(self, parent):
        tk_table_frame = ttk.Frame(parent)
        tk_table_frame.pack(fill="both", expand=True, pady=(5, 0))

        # column name in matches_table: (heading, width, anchor), this defines the order
        columns: dict[str, tuple[str, int, str]] = {
            "feature_id_meas": ("Meas. ID", 60, "center"),
            "feature_id_lib": ("Lib. ID", 60, "center"),
            "RT": ("RT (min)", 60, "center"),
            "name": ("Compound name", 180, "w"),
            "formula": ("Formula", 110, "w"),
            "ms2_score": ("MS2 score", 70, "e"),
            "n_hits_ms2": ("#MS2 peaks", 75, "e"),
            "dmz_mda": ("Δm/z (mDa)", 80, "e"),
            "dmz_ppm": ("Δm/z (ppm)", 80, "e"),
            "source_library": ("Source", 150, "w"),
        }
        self.table_columns: list[str] = list(columns)
        self.tk_matches_table: ttk.Treeview = ttk.Treeview(
            tk_table_frame,
            columns=self.table_columns,
            show="headings",
            selectmode="browse"
        )
        for column, (heading, width, anchor) in columns.items():
            self.tk_matches_table.heading(column, text=heading)
            self.tk_matches_table.column(column, width=width, anchor=anchor, stretch=column in ("name", "source_library"))
        self.tk_matches_table.bind("<<TreeviewSelect>>", self.on_select_table_row)

        y_scroll = ttk.Scrollbar(tk_table_frame, orient="vertical", command=self.tk_matches_table.yview)
        x_scroll = ttk.Scrollbar(tk_table_frame, orient="horizontal", command=self.tk_matches_table.xview)
        self.tk_matches_table.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.tk_matches_table.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")

        tk_table_frame.rowconfigure(0, weight=1)
        tk_table_frame.columnconfigure(0, weight=1)

    def _create_match_plot_panel(self, parent):
        self.tk_match_plot_frame = ttk.LabelFrame(parent, text="Match overview (select a row in the table)")
        parent.add(self.tk_match_plot_frame, weight=3)

    def _create_scatter_plot_panel(self, parent):
        tk_scatter_panel = ttk.LabelFrame(parent, text="Scatter plot m/z vs RT (click a feature to filter the table)")
        parent.add(tk_scatter_panel, weight=2)

        tk_scatter_buttons = ttk.Frame(tk_scatter_panel, padding=5)
        tk_scatter_buttons.pack(fill="x")
        ttk.Checkbutton(tk_scatter_buttons, text="Show annotations", variable=self.tk_show_annotations,
                        command=self.redraw_scatter_plot).pack(side="left", padx=3)
        ttk.Button(tk_scatter_buttons, text="Highlight NH3 loss",
                   command=self.on_click_highlight_nh3_loss).pack(side="left", padx=3)
        ttk.Button(tk_scatter_buttons, text="Remove highlight",
                   command=self.on_click_remove_highlight).pack(side="left", padx=3)

        self.tk_scatter_plot_frame = ttk.Frame(tk_scatter_panel)
        self.tk_scatter_plot_frame.pack(fill="both", expand=True)

    # ------------------------------------------------------------------------------------------------------------------
    # Loading databases
    # ------------------------------------------------------------------------------------------------------------------
    def on_browse_meas_db(self):
        """Connect to a measured database, reset all results and start loading its MS2 spectra."""
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
        self._destroy_match_plot()

        self.neutral_loss_f_ids = None
        self.draw_scatter_plot()

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

    # ------------------------------------------------------------------------------------------------------------------
    # Matching
    # ------------------------------------------------------------------------------------------------------------------
    def on_click_find_matches(self):
        """Validate the matching parameters, search the library for all measured features and fill the table."""
        if self.meas_db is None or self.lib_db is None:
            messagebox.showerror("Error", "Please select both databases.")
            return

        require_ms2: bool = self.tk_require_ms2.get()
        if require_ms2:
            if self.ms2_loading:
                messagebox.showinfo("MS2 spectra", "MS2 spectra are still loading. Try again shortly.")
                return
            if self.ms2_loading_error or self.ms2_spectra is None:
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

        self._clear_matches_table()
        self.tk_status_field.config(text="Searching the library ...")
        self.tk_find_matches_button.config(state="disabled")
        self.tk_root.update_idletasks()
        try:
            self.matches = self.lib_db.find_matches(
                mzs=self.meas_db.mzs,
                max_dmz_ppm=max_dmz_ppm,
                max_dmz_da=max_dmz_da,
                max_ms2_dmz_da=max_ms2_dmz_da,
                metric=metric,
                min_ms2_score=min_ms2_score,
                return_nhits_ms2=True,
                ms2_spectra=self.ms2_spectra,
                require_ms2=require_ms2
            )
            self.matches_table = self._matches_to_table(self.matches)
        except Exception as error:
            self.matches = {}
            self.matches_table = pd.DataFrame(columns=self.table_columns)
            self.filtered_matches_table = self.matches_table.copy()
            messagebox.showerror("Error", f"Search failed:\n{error}")
            self.tk_status_field.config(text="Failed while looking for matches.")
            return
        finally:
            self.tk_find_matches_button.config(state="normal")

        self.apply_filters(show_empty_message=False)
        count: int = len(self.matches_table)
        self.tk_status_field.config(text=f"Found {count} matches.")
        if count == 0:
            messagebox.showinfo("Search complete", "No matches found.")

    def _matches_to_table(self, matches: dict[int, list[dict]]) -> pd.DataFrame:
        """Turn the output of Library.find_matches into the table shown in the GUI (rounded, RT in minutes)."""
        table: pd.DataFrame = lib_match_result_to_table(matches)
        table.loc[:, 'ms2_score'] = table['ms2_score'].round(3)
        table.loc[:, 'dmz_mda'] = table['dmz_mda'].round(1).abs()
        table.loc[:, 'dmz_ppm'] = table['dmz_ppm'].round(1).abs()
        table.loc[:, 'RT'] = table.feature_id_meas.apply(
            lambda f_id: self.meas_db.retention_times_in_seconds[f_id] / 60
        ).round(2)
        return table

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
            self.filtered_matches_table = self.matches_table.copy()
            self._show_matches_in_table(self.filtered_matches_table)
            return

        table: pd.DataFrame = self.matches_table
        mask = np.ones(table.shape[0], dtype=bool)

        # numeric filters: (entry, column or series to compare, type, is minimum, error message)
        mzs_meas = table.feature_id_meas.map(self.meas_db.mzs)
        numeric_filters = [
            (self.tk_filter_min_f_id_meas, table.feature_id_meas, int, True, 'Minimum ID must be a whole number.'),
            (self.tk_filter_max_f_id_meas, table.feature_id_meas, int, False, 'Maximum ID must be a whole number.'),
            (self.tk_filter_min_rt, table.RT, float, True, 'Minimum RT must be a number.'),
            (self.tk_filter_max_rt, table.RT, float, False, 'Maximum RT must be a number.'),
            (self.tk_filter_min_mz, mzs_meas, float, True, 'Minimum m/z must be a number.'),
            (self.tk_filter_max_mz, mzs_meas, float, False, 'Maximum m/z must be a number.'),
            (self.tk_filter_min_ms2_score, table.ms2_score, float, True, 'Minimum MS2 score must be a number.'),
            (self.tk_filter_min_n_hits_ms2, table.n_hits_ms2, int, True,
             'Minimum #MS2 peak matches must be a whole number.'),
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

        self.filtered_matches_table = table.loc[mask, :].reset_index(drop=True)
        self._show_matches_in_table(self.filtered_matches_table)
        self.mark_filtered_features()
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

    def on_click_remove_filters(self):
        """Empty all filter fields and show all matches."""
        for tk_entry in self.tk_filter_entries:
            tk_entry.delete(0, tk.END)
        self.last_filter_values = self._get_filter_values()
        self.filtered_matches_table = self.matches_table.copy()
        self._show_matches_in_table(self.filtered_matches_table)
        self.mark_filtered_features()
        self.tk_status_field.config(text=f"Showing all {len(self.matches_table)} matches.")

    # ------------------------------------------------------------------------------------------------------------------
    # Results table
    # ------------------------------------------------------------------------------------------------------------------
    def _clear_matches_table(self):
        self.tk_matches_table.delete(*self.tk_matches_table.get_children())

    def _show_matches_in_table(self, table: pd.DataFrame):
        """Replace the rows of the results table, rows of the same measured feature share a background color."""
        self._clear_matches_table()

        f_ids_meas = sorted(table['feature_id_meas'].unique())
        n_colors: int = len(self.table_row_colors)
        for color_idx, f_id_meas in enumerate(f_ids_meas):
            self.tk_matches_table.tag_configure(
                f"feature_{f_id_meas}", background=self.table_row_colors[color_idx % n_colors]
            )

        for _, row in table.iterrows():
            values: list = ["" if pd.isna(value := row.get(column, "")) else value for column in self.table_columns]
            self.tk_matches_table.insert("", 'end', values=values, tags=(f"feature_{row.feature_id_meas}",))

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
        self.request_match_plot(int(match.feature_id_meas), int(match.feature_id_lib))

    def on_click_annotate(self):
        """Store the selected library match as a manual annotation in the measured database."""
        if (match := self._get_selected_match()) is None:
            messagebox.showerror("Error", "Please select a match row first.")
            return

        f_id_meas, f_id_lib = int(match.feature_id_meas), int(match.feature_id_lib)
        if not messagebox.askyesno("Confirm", f"Assign library feature {f_id_lib} to measured feature {f_id_meas}?"):
            return
        try:
            self.meas_db.add_annotations_from_library(
                library=self.lib_db,
                f_id_meas=f_id_meas,
                f_id_lib=f_id_lib,
                ms2_score=match.ms2_score,
            )
        except Exception as error:
            messagebox.showerror("Annotation Error", str(error))
            return
        messagebox.showinfo("Annotation", f"Annotation added for measured feature {f_id_meas}.")

    # ------------------------------------------------------------------------------------------------------------------
    # Match overview
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def _pack_canvas(canvas: FigureCanvasTkAgg):
        """Fill the parent with the canvas, the figure is scaled to the available space."""
        tk_widget = canvas.get_tk_widget()
        # without this, the canvas requests the full figure size and pushes the dividers of the panes
        tk_widget.config(width=1, height=1)
        tk_widget.pack(fill=tk.BOTH, expand=True)

    def request_match_plot(self, f_id_meas: int, f_id_lib: int):
        """
        Plot the match overview of a measured and library feature.

        Building the figure is slow, so it is done in a worker thread. Only one worker runs at a time;
        if the selection changes while plotting, only the latest request is plotted afterwards.
        """
        match_result: dict = next(m for m in self.matches[f_id_meas] if m['feature_id'] == f_id_lib)
        try:
            max_dmz_ppm: float = self.tk_max_dmz_ppm.get()
            max_dmz_da: float = self.tk_max_dmz_mda.get() * 1e-3
        except tk.TclError:
            messagebox.showerror("Invalid Tolerance", "Check the precursor tolerance values.")
            return
        # mass accuracy for the theoretical isotope pattern
        mz_tol_da: float = max(max_dmz_da, max_dmz_ppm * 1e-6 * self.meas_db.mzs[f_id_meas])

        self.pending_match_plot = dict(
            f_id_lib=f_id_lib,
            f_id_meas=f_id_meas,
            meas=self.meas_db,
            mz_tol_da=mz_tol_da,
            annotation_relative_cutoff=0.3,  # TODO: turn this into slider
            match_result=match_result,
        )
        self.tk_status_field.config(text=f"Plotting match {f_id_meas} -> {f_id_lib} ...")
        if not self.match_plot_running:
            self._start_match_plot_worker()

    def _start_match_plot_worker(self):
        plot_kwargs, self.pending_match_plot = self.pending_match_plot, None
        self.match_plot_running = True
        threading.Thread(
            target=self._match_plot_worker,
            args=(self.lib_db, plot_kwargs),
            daemon=True
        ).start()

    def _match_plot_worker(self, lib_db: Library, plot_kwargs: dict):
        """Worker thread: build the match overview figure (not attached to tk yet)."""
        fig = Figure(figsize=(15, 15))
        try:
            lib_db.plot_match(fig=fig, **plot_kwargs)
            error_text = None
        except Exception as error:
            error_text = str(error)

        try:
            self.tk_root.after(0, lambda: self._finish_match_plot(fig, plot_kwargs, error_text))
        except (tk.TclError, RuntimeError):  # window was closed in the meantime
            pass

    def _finish_match_plot(self, fig: Figure, plot_kwargs: dict, error_text: str | None):
        """Main thread: show the figure built by the worker, or start the next worker if the selection changed."""
        self.match_plot_running = False
        # selection changed while plotting, skip this figure and plot the latest selection instead
        if self.pending_match_plot is not None:
            self._start_match_plot_worker()
            return

        if error_text is not None:
            self.tk_status_field.config(text="Could not plot the selected match.")
            messagebox.showerror("Plot Error", f"Could not plot the selected match:\n{error_text}")
            return

        self._destroy_match_plot()
        self.match_canvas = FigureCanvasTkAgg(fig, master=self.tk_match_plot_frame)
        self.match_toolbar = NavigationToolbar2Tk(self.match_canvas, self.tk_match_plot_frame)
        self.match_toolbar.update()
        self.match_toolbar.pack()
        self._pack_canvas(self.match_canvas)
        self.match_canvas.draw()
        self.tk_status_field.config(text=f"Plotted match {plot_kwargs['f_id_meas']} -> {plot_kwargs['f_id_lib']}.")

    def _destroy_match_plot(self):
        if self.match_canvas is not None:
            self.match_canvas.get_tk_widget().destroy()
            self.match_canvas = None
        if self.match_toolbar is not None:
            self.match_toolbar.destroy()
            self.match_toolbar = None

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

    def draw_scatter_plot(self):
        """(Re)draw the scatter plot of the measurement including annotations, highlights and markers."""
        self._destroy_scatter_plot()

        fig = Figure(figsize=(7, 6))
        self.scatter_ax = fig.add_subplot(111)
        try:
            self.meas_db.plot_scatter(
                add_annotations=self.tk_show_annotations.get(),
                highlight_features=self.neutral_loss_f_ids,
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

        self.mark_filtered_features()
        if (match := self._get_selected_match()) is not None:
            self.mark_selected_feature(int(match.feature_id_meas))
        self.scatter_canvas.draw()

    def redraw_scatter_plot(self):
        """Redraw the scatter plot, if there is one (i.e. a measured database is loaded)."""
        if self.scatter_ax is not None:
            self.draw_scatter_plot()

    def on_click_highlight_nh3_loss(self):
        """Highlight features whose MS2 spectrum contains a neutral loss of NH3."""
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        try:
            self.neutral_loss_f_ids = self.meas_db.filter_by_neutral_loss(
                neutral_loss=CompoundDict('NH3').mass, mz_tol_da=5e-3
            )
        except Exception as error:
            messagebox.showerror("Error", f"Could not apply neutral loss filter:\n{error}")
            return
        self.draw_scatter_plot()
        self.tk_status_field.config(text=f"Highlighted {len(self.neutral_loss_f_ids)} features with NH3 loss.")

    def on_click_remove_highlight(self):
        self.neutral_loss_f_ids = None
        self.redraw_scatter_plot()

    def mark_filtered_features(self):
        """Circle features of the filtered table in black (only if a filter removes rows)."""
        if self.scatter_ax is None:
            return
        for artist in (self.filtered_features_marker, self.selected_feature_marker):
            if artist is not None:
                artist.remove()
        self.filtered_features_marker = self.selected_feature_marker = None

        if 0 < len(self.filtered_matches_table) < len(self.matches_table):
            points = [self._get_rt_mz(f_id) for f_id in self.filtered_matches_table.feature_id_meas.unique()]
            points = [p for p in points if p is not None]
            if points:
                rts, mzs = zip(*points)
                self.filtered_features_marker = self.scatter_ax.scatter(
                    rts, mzs, s=120, facecolors='none', edgecolors='black', linewidths=2, zorder=10,
                    label='Library match'
                )
                self.scatter_ax.legend(loc='upper right')
        self.scatter_canvas.draw_idle()

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
        self.scatter_ax.legend(loc='upper right')
        self.scatter_canvas.draw_idle()

    def on_click_scatter(self, event):
        """Filter the table to the feature closest to the click (only highlighted features if there are any)."""
        # ignore clicks outside the axes or while zooming / panning with the toolbar
        if event.inaxes is not self.scatter_ax or self.scatter_toolbar.mode:
            return

        candidate_f_ids = self.neutral_loss_f_ids if self.neutral_loss_f_ids else self.meas_db.feature_ids
        points = {f_id: p for f_id in candidate_f_ids if (p := self._get_rt_mz(f_id)) is not None}
        if not points:
            return
        f_ids = np.array(list(points.keys()))
        rts, mzs = np.array(list(points.values())).T

        # RT and m/z have different scales, so measure distance relative to the visible axis range
        x_min, x_max = self.scatter_ax.get_xlim()
        y_min, y_max = self.scatter_ax.get_ylim()
        distances = ((rts - event.xdata) / (x_max - x_min)) ** 2 + ((mzs - event.ydata) / (y_max - y_min)) ** 2
        closest_f_id = int(f_ids[np.argmin(distances)])

        for tk_entry in (self.tk_filter_min_f_id_meas, self.tk_filter_max_f_id_meas):
            tk_entry.delete(0, tk.END)
            tk_entry.insert(0, str(closest_f_id))
        self.apply_filters(show_empty_message=False)
        self.mark_selected_feature(closest_f_id)

        rt, mz = points[closest_f_id]
        self.tk_status_field.config(
            text=f"Selected feature {closest_f_id} (RT {rt:.2f} min, m/z {mz:.4f}) "
                 f"with {len(self.filtered_matches_table)} matches."
        )


if __name__ == "__main__":
    root: tk.Tk = tk.Tk()
    app: AnnotationGUI = AnnotationGUI(root)
    app.tk_root.mainloop()
