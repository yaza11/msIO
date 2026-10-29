"""
GUI for annotating measured features with matches from a spectral library.

Workflow:
    1. Select a measured database (FeatureManagerDB) and one or more library databases (Library).
       The scatter plot (m/z vs RT) of the measurement is drawn right away, the MS2 spectra of the measurement and
       the feature properties shown in the table (incl. I1 / I0 from the MS1 spectra) are loaded in the background.
    2. Set the matching parameters and click "Find Matches".
    3. Narrow down the results table with the filters (applied while typing, number fields take intervals like
       400-500; 512.3±0.01; the m/z can be set from a compound or an IPL batch of LipidCalculator), limit the rows per feature and add a dummy row per feature (only the feature
       properties, also for features without matches). Clicking a
       feature in the scatter plot selects its best match in the table. In the scatter plot, features with specific
       MS2 fragments or neutral losses can be highlighted and pairs of features with
       specific mass differences (e.g. adduct exchanges) can be connected.
    4. The overview panel shows the selected match (table row), a single feature (click in
       the scatter plot or feature ID field) or two compared features (shift-click two
       features or enter two IDs): MS2 spectra (peaks can be annotated by clicking them),
       MS1 spectrum, intensities in the samples and the structure (see feature_overview).
    5. Collect annotations in a local list: library matches ("Add library match"),
       compounds entered by hand or built with the LipidCalculator IPL builder ("Add manual
       annotation") and relations to other features: in-source fragments, duplicates, isotope peaks and other
       adducts ("Mark fragment / duplicate / isotope / adduct"), which can be hidden from the table.
       "Manage annotations" opens a window to review new and existing annotations and to
       write the changes to the measured database. "Export annotations as library" writes the annotated features to
       a new database that can be loaded as library.

Each area of the window is a collapsible section in its own module (*_panel.py), AnnotationGUI holds the shared
state (databases, matches, annotations, compared features) and connects the sections.

Slow work (loading MS2 spectra, searching the library, building the overview figure) runs in worker
threads. Worker threads never touch tk widgets; they hand their results back to the
main thread with ``tk_root.after``.
"""

import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Callable

import pandas as pd

from msIO.feature_managers.db import FeatureManagerDB, Library
from msIO.gui.annotation.annotations import AnnotationStore
from msIO.gui.annotation.common import SQLITE_FILE_TYPES, AnnotationKey
from msIO.gui.annotation.database_panel import DatabasePanel
from msIO.gui.annotation.dialogs import (
    AnnotationManagerWindow, FeatureRelationDialog, ManualAnnotationDialog, MzFromCompoundDialog
)
from msIO.gui.annotation.export import EXPORTED_KINDS, write_annotations_as_library
from msIO.gui.annotation.filter_panel import FilterError, FilterPanel
from msIO.gui.annotation.matching import (
    MzTarget, build_feature_table, combine_with_dummies, dummy_table, empty_matches_table, limit_matches_per_feature,
    match_targets, matches_to_table
)
from msIO.gui.annotation.matching_panel import MatchingPanel
from msIO.gui.annotation.overview_panel import OverviewPanel
from msIO.gui.annotation.results_panel import ResultsPanel
from msIO.gui.annotation.scatter_panel import ScatterPanel


class AnnotationGUI:
    """
    Main window for matching measured features against a library and annotating them.

    Layout:
        left:   databases, matching parameters, filters, results table (with the annotation buttons)
        right:  overview above the scatter plot, separated by a movable divider
        bottom: status bar
    All sections can be collapsed by clicking their header.

    Naming conventions:
        tk_*        tk widgets and tk variables
        on_*        callbacks of widgets and events
        _*worker    functions running in a worker thread
    """

    # height of the overview as fraction of the plot area, until the divider is moved
    OVERVIEW_HEIGHT_FRACTION: float = 0.55

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

        self.meas_db: FeatureManagerDB = None
        # loaded libraries by file name (stored with library annotations, so it identifies the library)
        self.lib_dbs: dict[str, Library] = {}
        # incremented for every loaded measured database, so results of workers for an older database are dropped
        self.meas_db_load_number: int = 0

        # MS2 spectra and properties of the measured features (see matching.build_feature_table), loaded in a worker
        self.ms2_spectra: dict = None
        self.feature_table: pd.DataFrame = None
        self.ms2_loading: bool = False
        self.ms2_loading_error: str = None

        # matching results: raw output of Library.find_matches (each match with the file name of its library under
        # "library") and as table (the filtered one is in the results panel)
        self.matches: dict[int, list[dict]] = {}
        self.matches_table: pd.DataFrame = empty_matches_table(list(ResultsPanel.COLUMNS))
        # matches with a dummy row per feature, built when needed: ((matches table, feature table), table)
        self._table_with_dummies: tuple[tuple, pd.DataFrame] = None
        # compounds of the m/z filter (e.g. an IPL batch) and the text they put into the m/z filter field, they are
        # dropped when the field is changed
        self.mz_targets: list[MzTarget] = []
        self.mz_targets_text: str = ""

        # manual annotations: changes are collected locally and only written to the measured database on request
        self.annotations: AnnotationStore = AnnotationStore()
        self.annotation_window: AnnotationManagerWindow = None

        self.compared_f_ids: list[int] = []  # up to two features selected by shift-click
        self.last_clicked_f_id: int = None  # in the scatter plot

        # position of the divider between overview and scatter plot while both are expanded
        self.plot_sash_position: int = None
        self.plot_sections_expanded: tuple[bool, bool] = (True, True)  # (overview, scatter plot)

        self._create_gui()

    # ------------------------------------------------------------------------------------------------------------------
    # GUI layout
    # ------------------------------------------------------------------------------------------------------------------
    def _create_gui(self):
        # status bar is packed before the main area, so it keeps its space when the window shrinks
        tk_status_bar = ttk.Frame(self.tk_root, relief="sunken")
        tk_status_bar.pack(side="bottom", fill="x")
        self.tk_status_field = ttk.Label(
            tk_status_bar, text="Please load both databases to begin.", anchor="w", padding=(5, 2)
        )
        self.tk_status_field.pack(side="left", fill="x", expand=True)
        # only shown while the library is searched or indexed, see show_progress
        self.tk_progress = ttk.Progressbar(tk_status_bar, mode="determinate", maximum=100, length=250)

        self.tk_main_area = ttk.PanedWindow(self.tk_root, orient="horizontal")
        self.tk_main_area.pack(fill="both", expand=True, padx=10, pady=(10, 5))

        # the results take the space left by the other sections (also when those are collapsed)
        tk_results_pane = ttk.Frame(self.tk_main_area)
        self.tk_main_area.add(tk_results_pane, weight=1)
        self.database = DatabasePanel(tk_results_pane, self)
        self.matching = MatchingPanel(tk_results_pane, self)
        self.filters = FilterPanel(tk_results_pane, self)
        self.results = ResultsPanel(tk_results_pane, self)
        for section in (self.database, self.matching, self.filters):
            section.pack(fill="x", pady=(0, 4))
        self.results.pack(fill="both", expand=True)

        self.tk_plot_area = ttk.PanedWindow(self.tk_main_area, orient="vertical")
        self.tk_main_area.add(self.tk_plot_area, weight=1)
        self.overview = OverviewPanel(self.tk_plot_area, self)
        self.scatter = ScatterPanel(self.tk_plot_area, self)
        self.tk_plot_area.add(self.overview, weight=3)
        self.tk_plot_area.add(self.scatter, weight=2)

        # the panes would otherwise be sized by the requested size of their content, so place the dividers once the
        # window is shown: results | plots at half the width, overview | scatter plot at OVERVIEW_HEIGHT_FRACTION
        def set_initial_dividers(event):
            if event.width <= 1:
                return
            self.tk_main_area.unbind("<Configure>")
            self.tk_main_area.sashpos(0, event.width // 2)
            self.tk_plot_area.sashpos(0, int(event.height * self.OVERVIEW_HEIGHT_FRACTION))

        self.tk_main_area.bind("<Configure>", set_initial_dividers)

    def on_toggle_plot_section(self):
        """
        Overview or scatter plot collapsed / expanded: a collapsed section only keeps its header, the other one gets
        the space. Expanding both again restores the divider.
        """
        overview_expanded, scatter_expanded = self.overview.expanded, self.scatter.expanded
        # the divider hasn't moved yet: if both were expanded before this toggle, remember where it was
        if self.plot_sections_expanded == (True, True):
            self.plot_sash_position = self.tk_plot_area.sashpos(0)
        self.plot_sections_expanded = (overview_expanded, scatter_expanded)

        # a collapsed section must not grow when the window is resized
        self.tk_plot_area.pane(self.overview, weight=3 if overview_expanded else 0)
        self.tk_plot_area.pane(self.scatter, weight=2 if scatter_expanded else 0)
        self.tk_root.update_idletasks()
        height: int = self.tk_plot_area.winfo_height()
        if overview_expanded and scatter_expanded:
            position = self.plot_sash_position or int(height * self.OVERVIEW_HEIGHT_FRACTION)
        elif not overview_expanded:
            position = self.overview.header_height()
        else:
            position = height - self.scatter.header_height()
        self.tk_plot_area.sashpos(0, position)

    # ------------------------------------------------------------------------------------------------------------------
    # Status bar and threads
    # ------------------------------------------------------------------------------------------------------------------
    def set_status(self, text: str):
        self.tk_status_field.config(text=text)

    def call_in_main_thread(self, func: Callable[[], None]):
        """Hand a function from a worker thread to the tk main loop."""
        try:
            self.tk_root.after(0, func)
        except (tk.TclError, RuntimeError):  # window was closed in the meantime
            pass

    def show_progress(self, text: str, indeterminate: bool = False):
        """Show the progress bar in the status bar, indeterminate: moving bar for work without known progress."""
        self.set_status(text)
        self.tk_progress.config(value=0, mode="indeterminate" if indeterminate else "determinate")
        self.tk_progress.pack(side="right", padx=5, before=self.tk_status_field)
        if indeterminate:
            self.tk_progress.start(20)

    def update_progress(self, percent: float, text: str):
        self.tk_progress.config(value=percent)
        self.set_status(text)

    def hide_progress(self):
        self.tk_progress.stop()
        self.tk_progress.pack_forget()

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
        self.database.tk_meas_path.set(path)
        self.meas_db_load_number += 1

        # clear data belonging to the previous database
        self.ms2_spectra = None
        self.feature_table = None
        self.ms2_loading_error = None
        self.ms2_loading = True
        self.matches = {}
        self.matches_table = empty_matches_table(self.results.columns)
        self.results.clear()
        self.annotations.reset(self.meas_db)
        self.compared_f_ids = []
        self.last_clicked_f_id = None
        self.overview.clear()
        self.overview.peak_annotations = {}
        self.scatter.reset()
        self._on_annotations_changed()

        self.set_status(f"Measured database loaded, loading spectra ...: {path}")
        threading.Thread(
            target=self._ms2_spectra_worker,
            args=(self.meas_db, self.meas_db_load_number),
            daemon=True
        ).start()

    def _ms2_spectra_worker(self, meas_db: FeatureManagerDB, load_number: int):
        """Worker thread: load all MS2 spectra of the measurement and the properties of the features."""
        spectra, feature_table, errors = None, None, []
        try:
            spectra = meas_db.get_ms_spectra(meas_db.feature_ids, level=2)
        except Exception as error:
            errors.append(f"MS2 spectra: {error}")
        try:
            # the MS1 spectra are only needed for I1 / I0
            feature_table = build_feature_table(meas_db, meas_db.get_ms_spectra(meas_db.feature_ids, level=1))
        except Exception as error:
            errors.append(f"feature properties: {type(error).__name__}: {error}")
        error_text: str | None = "\n".join(errors) or None
        self.call_in_main_thread(lambda: self._finish_ms2_loading(load_number, spectra, feature_table, error_text))

    def _finish_ms2_loading(self, load_number: int, spectra: dict | None, feature_table: pd.DataFrame | None,
                            error_text: str | None):
        """Main thread: store the loaded data, unless another database was loaded in the meantime."""
        if load_number != self.meas_db_load_number:
            return
        self.ms2_loading = False
        self.ms2_spectra = spectra
        self.feature_table = feature_table
        self.ms2_loading_error = None if spectra is not None else error_text

        if error_text is not None:
            self.set_status("Measured database loaded, but not all data could be loaded.")
            messagebox.showerror("Loading error", f"Could not load:\n{error_text}")
        else:
            self.set_status("Measured database, spectra and feature properties loaded.")
        if self.filters.tk_add_dummies.get():
            self.apply_filters(show_empty_message=False)

    def on_browse_lib_db(self):
        """Replace the libraries by the selected ones."""
        if paths := filedialog.askopenfilenames(filetypes=SQLITE_FILE_TYPES):
            self.load_libraries(list(paths))

    def on_add_lib_db(self):
        """Add the selected libraries to the loaded ones."""
        if paths := filedialog.askopenfilenames(filetypes=SQLITE_FILE_TYPES):
            self.load_libraries([lib_db.path_file for lib_db in self.lib_dbs.values()] + list(paths))

    def load_libraries(self, paths: list[str]):
        """
        Use the libraries at paths (libraries that are already loaded are kept). Library annotations are identified by
        the file name of the library, so two libraries can't have the same file name.
        """
        loaded: dict[str, Library] = {os.path.normcase(os.path.abspath(lib_db.path_file)): lib_db
                                      for lib_db in self.lib_dbs.values()}
        lib_dbs: dict[str, Library] = {}
        for path in dict.fromkeys(paths):  # without duplicates, in order
            name: str = os.path.basename(path)
            if name in lib_dbs:
                messagebox.showerror("Error", f"Two libraries are called {name}, please rename one of them.")
                return
            try:
                lib_dbs[name] = loaded.get(os.path.normcase(os.path.abspath(path))) or Library(path)
            except Exception as error:
                messagebox.showerror("Error", f"Failed to load library database {path}:\n{error}")
                return
        changed: bool = list(lib_dbs.values()) != list(self.lib_dbs.values())
        new: list[Library] = [lib_db for lib_db in lib_dbs.values() if lib_db not in self.lib_dbs.values()]
        self.lib_dbs = lib_dbs
        self.database.set_lib_paths([lib_db.path_file for lib_db in lib_dbs.values()])
        if not changed:
            return
        self.set_status(f"{len(lib_dbs)} library database(s) loaded: {', '.join(lib_dbs)}")
        self._clear_results_of_previous_library()
        self._offer_missing_indexes(new)

    def get_match_result(self, f_id_meas: int, library_name: str, f_id_lib: int) -> dict:
        """Match of a measured feature with a library feature, as returned by Library.find_matches."""
        return next(match for match in self.matches[f_id_meas]
                    if match['library'] == library_name and match['feature_id'] == f_id_lib)

    def _clear_results_of_previous_library(self):
        """Matches refer to the previous libraries, which may not be loaded anymore: search again."""
        if self.meas_db is None:
            return
        self.matches = {}
        self.matches_table = empty_matches_table(self.results.columns)
        self.apply_filters(show_empty_message=False)  # dummy rows stay
        if self.overview.shown_key is not None and self.overview.shown_key[0] == "match":
            self.overview.clear()  # also drops a match overview that is still being built

    def _offer_missing_indexes(self, dbs: list[FeatureManagerDB]):
        """
        Databases created before the indexes were added to the models load spectra much slower (full-scale library:
        search 66 s instead of 14 s). Creating them changes the files, so ask first.
        """
        without_index: list[FeatureManagerDB] = []
        for db in dbs:
            try:
                if db.get_missing_indexes():
                    without_index.append(db)
            except Exception:  # e.g. not an msIO database, the search will report it
                pass
        if not without_index:
            return
        sizes: str = ", ".join(f"{os.path.basename(db.path_file)} {os.path.getsize(db.path_file) / 1e9:.1f} GB"
                               for db in without_index)
        if not messagebox.askyesno(
                "Speed up the library",
                f"{len(without_index)} library database(s) have no index for loading spectra, which makes the search "
                f"several times slower.\n\n"
                f"Add the index now? This is done once and changes the file: it takes a while for large libraries "
                f"(about 20 s per 5 GB) and makes the file about 20 % larger (currently {sizes})."
        ):
            return
        self.matching.set_searching(True)
        self.show_progress("Creating the index of the library ...", indeterminate=True)

        def worker():
            try:
                for db in without_index:
                    db.create_missing_indexes()
                error_text = None
            except Exception as error:
                error_text = f"{type(error).__name__}: {error}"
            self.call_in_main_thread(lambda: finish(error_text))

        def finish(error_text: str | None):
            self.matching.set_searching(False)
            self.hide_progress()
            if error_text is not None:
                messagebox.showerror("Index", f"Could not create the index:\n{error_text}")
                self.set_status("Library loaded without index.")
            else:
                self.set_status("Index of the library created.")

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
        if self.meas_db is None or not self.lib_dbs:
            messagebox.showerror("Error", "Please select a measured database and at least one library.")
            return
        try:
            parameters: dict = self.matching.get_search_parameters()
        except ValueError as error:
            messagebox.showerror("Invalid input", str(error))
            return
        # without spectra, all MS2 scores would be missing
        if self.ms2_loading:
            messagebox.showinfo("MS2 spectra", "MS2 spectra are still loading. Try again shortly.")
            return
        if self.feature_table is None:
            messagebox.showerror("Error", "The properties of the measured features could not be loaded.")
            return
        if parameters['require_ms2'] and (self.ms2_loading_error or self.ms2_spectra is None):
            messagebox.showerror("MS2 spectra", "MS2 spectra could not be loaded.")
            return

        self.matching.set_searching(True)
        self.show_progress("Searching the library ...")
        search_kwargs: dict = dict(
            mzs=self.meas_db.mzs,
            return_nhits_ms2=True,
            # empty dict: scores are missing, but find_matches still returns the matches as dicts
            ms2_spectra=self.ms2_spectra or {},
            **parameters,
        )
        threading.Thread(
            target=self._find_matches_worker,
            args=(self.feature_table, dict(self.lib_dbs), search_kwargs, self.meas_db_load_number),
            daemon=True
        ).start()

    def _find_matches_worker(self, feature_table: pd.DataFrame, lib_dbs: dict[str, Library], search_kwargs: dict,
                             load_number: int):
        """Worker thread: search the libraries one after another and build the results table."""
        last_percent: list[int] = [-1]

        def report_progress(idx_library: int, n_done: int, n_total: int):
            # only hand over changes of whole percents, the main thread doesn't need every feature
            percent = int((idx_library + n_done / max(n_total, 1)) / len(lib_dbs) * 100)
            if percent != last_percent[0]:
                last_percent[0] = percent
                text = f"Scoring MS2 matches (library {idx_library + 1} of {len(lib_dbs)}) ... {percent} %"
                self.call_in_main_thread(lambda: self.update_progress(percent, text))

        matches, table, error_text = None, None, None
        try:
            matches = {}
            for idx_library, (library_name, lib_db) in enumerate(lib_dbs.items()):
                found = lib_db.find_matches(**search_kwargs, progress_callback=lambda n_done, n_total, i=idx_library:
                                            report_progress(i, n_done, n_total))
                for f_id_meas, matches_of_feature in found.items():
                    for match in matches_of_feature:
                        match['library'] = library_name
                    matches.setdefault(f_id_meas, []).extend(matches_of_feature)
            self.call_in_main_thread(lambda: self.update_progress(100, "Building the results table ..."))
            table = matches_to_table(matches, feature_table, lib_dbs, search_kwargs['max_dmz_ppm'],
                                     search_kwargs['max_dmz_da'], list(ResultsPanel.COLUMNS))
        except Exception as error:
            error_text = f"{type(error).__name__}: {error}"
        self.call_in_main_thread(lambda: self._finish_find_matches(load_number, lib_dbs, matches, table, error_text))

    def _finish_find_matches(self, load_number: int, lib_dbs: dict[str, Library], matches: dict | None,
                             table: pd.DataFrame | None, error_text: str | None):
        """Main thread: show the results, unless another database was loaded during the search."""
        self.matching.set_searching(False)
        self.hide_progress()
        if load_number != self.meas_db_load_number or lib_dbs != self.lib_dbs:
            self.set_status("Search results dropped, a different database was loaded meanwhile.")
            return
        if error_text is not None:
            messagebox.showerror("Error", f"Search failed:\n{error_text}")
            self.set_status("Failed while looking for matches.")
            return

        self.matches, self.matches_table = matches, table
        self.overview.shown_key = None  # scores may have changed, selecting the match again replots it
        self.apply_filters(show_empty_message=False)
        count: int = len(self.matches_table)
        self.set_status(f"Found {count} matches in {len(lib_dbs)} librar{'y' if len(lib_dbs) == 1 else 'ies'}.")
        if count == 0:
            messagebox.showinfo("Search complete", "No matches found.")

    # ------------------------------------------------------------------------------------------------------------------
    # Filters
    # ------------------------------------------------------------------------------------------------------------------
    def get_unfiltered_table(self) -> pd.DataFrame:
        """All rows of the results table: the matches and, if selected, a dummy row per measured feature."""
        if not self.filters.tk_add_dummies.get() or self.feature_table is None:
            return self.matches_table
        sources: tuple = (self.matches_table, self.feature_table)
        if self._table_with_dummies is None or any(a is not b for a, b in zip(self._table_with_dummies[0], sources)):
            dummies = dummy_table(self.feature_table, list(ResultsPanel.COLUMNS))
            self._table_with_dummies = (sources, combine_with_dummies(self.matches_table, dummies))
        return self._table_with_dummies[1]

    def apply_filters(self, show_empty_message: bool = True, show_dialogs: bool = True):
        """
        Filter the results (matches and dummy rows) by all filter fields, limit the rows per feature and show the
        result in the table.

        Invalid filter values leave the table unchanged. They are reported in a dialog,
        or in the status bar if show_dialogs is False (used while typing).
        """
        table: pd.DataFrame = self.get_unfiltered_table()
        if table.empty:
            self._show_filtered_matches(table.copy())
            return
        try:
            mask = self.filters.compute_mask(table)
            max_rows: int | None = self.filters.get_max_rows()
        except FilterError as error:
            # while typing, half-entered values are expected, so don't open a dialog for every key stroke
            if show_dialogs:
                messagebox.showerror(error.title, error.text)
            else:
                self.set_status(f"Filters not updated: {error.text}")
            return

        filtered: pd.DataFrame = table.loc[mask, :]
        if max_rows is not None:
            filtered = limit_matches_per_feature(filtered, max_rows)
        if self.mz_targets and self.filters.tk_mz.get().strip() != self.mz_targets_text:
            self.mz_targets, self.mz_targets_text = [], ""
        if self.mz_targets:
            filtered = filtered.assign(target=["; ".join(target.label for target in targets)
                                               for targets in self.get_row_targets(filtered.mz_meas)])
        self._show_filtered_matches(filtered)
        n_shown: int = len(self.results.table)
        if show_empty_message and n_shown == 0:
            self.set_status("No rows meet the current filters.")
        elif n_shown:
            self.set_status(f"Showing {n_shown} of {len(table)} rows.")

    def _show_filtered_matches(self, table: pd.DataFrame):
        self.results.set_table(table)
        self.scatter.mark_filtered_features()

    def get_filtered_feature_ids(self) -> list[int]:
        """Measured features of the shown rows, empty if the filters don't remove any row."""
        if 0 < len(self.results.table) < len(self.get_unfiltered_table()):
            return [int(f_id) for f_id in self.results.table.feature_id_meas.unique()]
        return []

    def get_row_targets(self, mzs) -> list[list[MzTarget]]:
        """Targets of the m/z filter at each m/z (precursor tolerance), closest first."""
        try:
            max_dmz_ppm, max_dmz_da = self.matching.precursor_tolerances()
        except tk.TclError:
            return [[] for _ in mzs]
        return match_targets(pd.to_numeric(pd.Series(mzs), errors='coerce').to_numpy(dtype=float), self.mz_targets,
                             max_dmz_ppm, max_dmz_da)

    def set_mz_targets(self, targets: list[MzTarget]):
        """
        Filter the results to the m/z of the targets (each ± precursor tolerance, e.g. an IPL batch), show the dummy
        rows (features without library matches) and name the fitting targets in the table.
        """
        try:
            mzs = sorted({round(target.mz, 4) for target in targets})
            text: str = "; ".join(f"{mz:.4f}±{self.matching.precursor_tolerance_da(mz):.4f}" for mz in mzs)
        except tk.TclError:
            messagebox.showerror("Invalid input", "Check the precursor tolerances.")
            return
        self.mz_targets, self.mz_targets_text = list(targets), text
        self.filters.set_mz_text(text)
        self.filters.tk_add_dummies.set(True)
        self.apply_filters(show_empty_message=False)
        n_features: int = len(set(self.results.table.feature_id_meas))
        self.set_status(f"m/z filter of {len(targets)} targets ({len(mzs)} m/z values): {n_features} features "
                        f"(dummy rows show the features without library matches, see column Target).")

    def set_mz_filter(self, mz: float):
        """Filter the results to measured features at mz (± precursor tolerance)."""
        try:
            tolerance: float = self.matching.precursor_tolerance_da(mz)
        except tk.TclError:
            messagebox.showerror("Invalid input", "Check the precursor tolerances.")
            return
        self.filters.set_mz_text(f"{mz:.4f}±{tolerance:.4f}")
        self.apply_filters()

    def on_click_mz_from_compound(self):
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        MzFromCompoundDialog(self)

    def on_click_remove_filters(self):
        """Empty all filter fields and show all rows."""
        self.filters.clear()
        self.apply_filters(show_empty_message=False)
        self.set_status(f"Showing all {len(self.results.table)} rows.")

    # ------------------------------------------------------------------------------------------------------------------
    # Manual annotations
    # ------------------------------------------------------------------------------------------------------------------
    def _confirm_discarding_unsaved_annotations(self) -> bool:
        to_add, to_remove = self.annotations.unsaved_changes()
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
        added: bool = self.annotations.add(key, info)
        self._on_annotations_changed()
        return added

    def _get_current_feature(self) -> int | None:
        """Feature to prefill annotation dialogs: selected in the table, else last clicked in the scatter plot."""
        if (match := self.results.selected_match()) is not None:
            return int(match.feature_id_meas)
        if self.last_clicked_f_id is not None:
            return self.last_clicked_f_id
        return self.compared_f_ids[-1] if self.compared_f_ids else None

    def on_click_add_library_match(self):
        """Add the selected library match to the local list of annotations."""
        if (match := self.results.selected_match()) is None or match.is_dummy:
            messagebox.showerror("Error", "Please select a row with a library match first.")
            return
        f_id_meas, f_id_lib, library_name = int(match.feature_id_meas), int(match.feature_id_lib), match.library
        match_result: dict = self.get_match_result(f_id_meas, library_name, f_id_lib)
        key: AnnotationKey = AnnotationStore.library_key(f_id_meas, library_name, f_id_lib)
        added: bool = self.add_pending_annotation(key, dict(
            kind="library",
            name=match_result.get('name'),
            formula=match_result.get('formula'),
            adduct=match_result.get('adduct'),
            ms2_score=match_result.get('ms2_score'),
            reference=f_id_lib,
            library_name=library_name,
            library=self.lib_dbs[library_name],
        ))
        if added:
            self.set_status(f"Added library feature {f_id_lib} ({library_name}) for feature {f_id_meas} to the "
                            f"annotations. Use \"Manage annotations\" to write them to the database.")
        else:
            self.set_status(f"Library feature {f_id_lib} ({library_name}) is already annotated for feature "
                            f"{f_id_meas}.")

    def on_click_add_manual_annotation(self):
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        # a target of the m/z filter at the selected row is the suggested compound
        compound: MzTarget | None = None
        if self.mz_targets and (match := self.results.selected_match()) is not None:
            compound = next(iter(self.get_row_targets([match.mz_meas])[0]), None)
        ManualAnnotationDialog(self, self._get_current_feature(), compound=compound)

    def on_click_mark_relation(self, kind: str = "in_source_fragment"):
        """
        Open the dialog to mark a feature as in-source fragment, duplicate, isotope peak or adduct of another one,
        prefilled with the compared features (the dialog orders them by m/z for fragments and isotope peaks).
        """
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        if len(self.compared_f_ids) == 2:
            f_id, f_id_reference = sorted(self.compared_f_ids, key=lambda f_id_: self.meas_db.mzs[f_id_])
        else:
            f_id, f_id_reference = self._get_current_feature(), None
        return FeatureRelationDialog(self, f_id, f_id_reference, kind=kind)

    def toggle_annotations(self, keys: list[AnnotationKey], keep: bool = None):
        """Mark annotations to keep or discard (keep=None toggles each of them)."""
        self.annotations.toggle(keys, keep)
        self._on_annotations_changed()

    def write_annotations_to_db(self, parent: tk.Misc = None):
        """Add the new and remove the discarded annotations in the measured database (after confirmation)."""
        to_add, to_remove = self.annotations.unsaved_changes()
        if not (to_add or to_remove):
            messagebox.showinfo("Annotations", "There are no changes to write.", parent=parent)
            return
        if not messagebox.askyesno(
                "Write annotations",
                f"Add {len(to_add)} and remove {len(to_remove)} annotation(s) in the measured database?",
                parent=parent
        ):
            return

        n_added, n_removed, errors = self.annotations.write(self.meas_db)
        self._on_annotations_changed()
        if self.scatter.tk_show_annotations.get():
            self.scatter.redraw()
        if errors:
            messagebox.showerror("Annotation Error", "\n".join(errors), parent=parent)
        self.set_status(f"Wrote annotations to the database: {n_added} added, {n_removed} removed.")

    def _on_annotations_changed(self):
        """Update everything that shows annotation states: table column, button text, manager window, scatter plot."""
        if self.filters.tk_hide_related.get():  # the hidden features changed
            self.apply_filters(show_empty_message=False, show_dialogs=False)
        self.results.refresh_annotation_column()
        self.results.set_unsaved_annotation_count(sum(len(keys) for keys in self.annotations.unsaved_changes()))
        if self.annotation_window is not None:
            self.annotation_window.refresh()
        self.scatter.mark_annotated_features()

    def on_click_export_library(self):
        """
        Write the annotated features (library matches and manual annotations as listed in "Manage annotations", incl.
        the unsaved changes) to a new database that can be loaded as library.
        """
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        annotations: dict[AnnotationKey, dict] = self.annotations.kept()
        if not (n_entries := sum(key[0] in EXPORTED_KINDS for key in annotations)):
            messagebox.showinfo("Export", "There are no library matches or manual annotations to export.")
            return
        stem: str = os.path.splitext(os.path.basename(self.meas_db.path_file))[0]
        path: str = filedialog.asksaveasfilename(
            title="Export annotations as library", defaultextension=".sqlite", filetypes=SQLITE_FILE_TYPES,
            initialfile=f"{stem}_annotations.sqlite"
        )
        if not path:
            return
        # the file is replaced, so it must not be one of the open databases
        open_paths = {os.path.normcase(os.path.abspath(db.path_file)) for db in [self.meas_db, *self.lib_dbs.values()]}
        if os.path.normcase(os.path.abspath(path)) in open_paths:
            messagebox.showerror("Export", "Please choose a new file, this database is in use.")
            return
        source_name: str = os.path.splitext(os.path.basename(path))[0]
        self.set_status(f"Exporting {n_entries} annotation(s) ...")
        self.tk_root.update_idletasks()
        try:
            n_written: int = write_annotations_as_library(path, self.meas_db, annotations, self.lib_dbs, source_name)
        except Exception as error:
            messagebox.showerror("Export", f"Could not export the annotations:\n{type(error).__name__}: {error}")
            self.set_status("Export failed.")
            return
        n_skipped: int = len(annotations) - n_written
        self.set_status(f"Exported {n_written} annotation(s) as library to {path}"
                        + (f" ({n_skipped} in-source fragment annotation(s) left out)." if n_skipped else "."))

    def on_click_manage_annotations(self):
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        if self.annotation_window is None:
            self.annotation_window = AnnotationManagerWindow(self)
        else:
            self.annotation_window.lift()

    # ------------------------------------------------------------------------------------------------------------------
    # Selected and compared features
    # ------------------------------------------------------------------------------------------------------------------
    def get_rt_mz(self, f_id: int) -> tuple[float, float] | None:
        """RT in minutes (as in FeatureManagerDB.plot_scatter) and m/z of a measured feature."""
        rt = self.meas_db.retention_times_in_seconds.get(f_id)
        mz = self.meas_db.mzs.get(f_id)
        if rt is None or mz is None:
            return None
        return rt / 60, mz

    def on_match_selected(self, match: pd.Series, by_program: bool = False):
        """Row of the results table selected: mark its feature in the scatter plot and show the match (or feature)."""
        self.scatter.mark_selected_feature(int(match.feature_id_meas))
        if by_program:  # selected by the program (scatter click, restored selection): keep the overview
            return
        # e.g. selected again after sorting or filtering
        if self._overview_key_of_row(match) != self.overview.shown_key:
            self.clear_comparison()
            self._show_row(match)

    @staticmethod
    def _overview_key_of_row(match: pd.Series) -> tuple:
        """Overview of a row of the results table (see OverviewPanel.shown_key): the match, the feature for dummies."""
        if match.is_dummy:
            return "feature", (int(match.feature_id_meas),), None
        return "match", (int(match.feature_id_meas),), (match.library, int(match.feature_id_lib))

    def _show_row(self, match: pd.Series):
        mode, f_ids, lib_ref = self._overview_key_of_row(match)
        self.overview.request(mode, list(f_ids), lib_ref)

    def on_scatter_click(self, f_id: int, shift: bool = False):
        """
        Click: select the best match of the feature in the results table and show it, features without
        matches (in the filtered table) are shown alone.
        Shift-click: add the feature to the comparison (the oldest of two is replaced).
        """
        self.last_clicked_f_id = f_id
        if shift:
            self.add_feature_to_comparison(f_id)
            return
        self.clear_comparison()
        # the overview is requested here, the selection event isn't generated if the row was already selected
        if self.results.select_best_match(f_id, show_match=False):
            self.scatter.mark_selected_feature(f_id)
            self._show_row(self.results.selected_match())
        else:
            self.show_feature(f_id)

    def show_feature(self, f_id: int):
        """Overview of a single measured feature (with or without matches)."""
        self.scatter.mark_selected_feature(f_id)
        self.overview.request("feature", [f_id])

    def show_feature_ids(self, f_ids: list[int]):
        """Feature ID field of the overview: one feature is shown alone, two are compared."""
        if len(f_ids) == 1:
            self.clear_comparison()
            self.show_feature(f_ids[0])
        else:
            self.compared_f_ids = list(f_ids)
            self.scatter.mark_compared_features()
            self.overview.request("comparison", list(f_ids))

    def clear_comparison(self):
        """The overview shows something else than the shift-clicked features: remove them and their squares."""
        if self.compared_f_ids:
            self.compared_f_ids = []
            self.scatter.mark_compared_features()

    def add_feature_to_comparison(self, f_id: int):
        """Add a feature to the comparison: the first one is shown alone, two are compared."""
        if f_id in self.compared_f_ids:
            self.set_status(f"Feature {f_id} is already compared.")
            return
        self.compared_f_ids = (self.compared_f_ids + [f_id])[-2:]
        self.scatter.mark_compared_features()
        if len(self.compared_f_ids) < 2:
            self.overview.request("feature", [f_id])
        else:
            self.overview.request("comparison", list(self.compared_f_ids))
