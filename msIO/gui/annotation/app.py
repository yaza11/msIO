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
    AnnotationManagerWindow, InSourceFragmentDialog, ManualAnnotationDialog, MzFromCompoundDialog
)
from msIO.gui.annotation.filter_panel import FilterError, FilterPanel
from msIO.gui.annotation.matching import empty_matches_table, matches_to_table
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
        self.lib_db: Library = None
        # incremented for every loaded measured database, so results of workers for an older database are dropped
        self.meas_db_load_number: int = 0

        # MS2 spectra of the measurement, loaded in a worker thread
        self.ms2_spectra: dict = None
        self.ms2_loading: bool = False
        self.ms2_loading_error: str = None

        # matching results: raw output of Library.find_matches and as table (the filtered one is in the results panel)
        self.matches: dict[int, list[dict]] = {}
        self.matches_table: pd.DataFrame = empty_matches_table(list(ResultsPanel.COLUMNS))

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

        self.set_status(f"Measured database loaded, loading MS2 spectra ...: {path}")
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
        self.call_in_main_thread(lambda: self._finish_ms2_loading(load_number, spectra, error_text))

    def _finish_ms2_loading(self, load_number: int, spectra: dict | None, error_text: str | None):
        """Main thread: store the loaded MS2 spectra, unless another database was loaded in the meantime."""
        if load_number != self.meas_db_load_number:
            return
        self.ms2_loading = False
        self.ms2_spectra = spectra
        self.ms2_loading_error = error_text

        if error_text is not None:
            self.set_status("Measured database loaded, but MS2 spectra could not be loaded.")
            messagebox.showerror("MS2 loading error", f"Could not load MS2 spectra:\n{error_text}")
            return
        self.set_status("Measured database and MS2 spectra loaded.")

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
        self.database.tk_lib_path.set(path)
        self.set_status(f"Library database loaded: {path}")
        self._clear_results_of_previous_library()
        self._offer_missing_indexes(self.lib_db)

    def _clear_results_of_previous_library(self):
        """Matches refer to library feature ids of the previous library, so they can't be shown with the new one."""
        if self.meas_db is None:
            return
        self.matches = {}
        self.matches_table = empty_matches_table(self.results.columns)
        self.results.clear()
        self.scatter.mark_filtered_features()
        if self.overview.shown_key is not None and self.overview.shown_key[0] == "match":
            self.overview.clear()  # also drops a match overview that is still being built

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
        self.matching.set_searching(True)
        self.show_progress("Creating the index of the library ...", indeterminate=True)

        def worker():
            try:
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
        if self.meas_db is None or self.lib_db is None:
            messagebox.showerror("Error", "Please select both databases.")
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
                self.call_in_main_thread(lambda: self.update_progress(percent, f"Scoring MS2 matches ... {percent} %"))

        matches, table, error_text = None, None, None
        try:
            matches = lib_db.find_matches(**search_kwargs, progress_callback=report_progress)
            self.call_in_main_thread(lambda: self.update_progress(100, "Building the results table ..."))
            table = matches_to_table(matches, meas_db, lib_db, search_kwargs['max_dmz_ppm'],
                                     search_kwargs['max_dmz_da'], list(ResultsPanel.COLUMNS))
        except Exception as error:
            error_text = f"{type(error).__name__}: {error}"
        self.call_in_main_thread(lambda: self._finish_find_matches(load_number, lib_db, matches, table, error_text))

    def _finish_find_matches(self, load_number: int, lib_db: Library, matches: dict | None, table: pd.DataFrame | None,
                             error_text: str | None):
        """Main thread: show the results, unless another database was loaded during the search."""
        self.matching.set_searching(False)
        self.hide_progress()
        if load_number != self.meas_db_load_number or lib_db is not self.lib_db:
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
        self.set_status(f"Found {count} matches.")
        if count == 0:
            messagebox.showinfo("Search complete", "No matches found.")

    # ------------------------------------------------------------------------------------------------------------------
    # Filters
    # ------------------------------------------------------------------------------------------------------------------
    def apply_filters(self, show_empty_message: bool = True, show_dialogs: bool = True):
        """
        Filter matches_table by all filter fields and show the result in the table.

        Invalid filter values leave the table unchanged. They are reported in a dialog,
        or in the status bar if show_dialogs is False (used while typing).
        """
        if self.matches_table.empty:
            self._show_filtered_matches(self.matches_table.copy())
            return
        try:
            mask = self.filters.compute_mask(self.matches_table)
        except FilterError as error:
            # while typing, half-entered values are expected, so don't open a dialog for every key stroke
            if show_dialogs:
                messagebox.showerror(error.title, error.text)
            else:
                self.set_status(f"Filters not updated: {error.text}")
            return

        self._show_filtered_matches(self.matches_table.loc[mask, :])
        n_shown: int = len(self.results.table)
        if show_empty_message and n_shown == 0:
            self.set_status("No matches meet the current filters.")
        elif n_shown:
            self.set_status(f"Showing {n_shown} of {len(self.matches_table)} matches.")

    def _show_filtered_matches(self, table: pd.DataFrame):
        self.results.set_table(table)
        self.scatter.mark_filtered_features()

    def get_filtered_feature_ids(self) -> list[int]:
        """Measured features of the shown matches, empty if the filters don't remove any match."""
        if 0 < len(self.results.table) < len(self.matches_table):
            return [int(f_id) for f_id in self.results.table.feature_id_meas.unique()]
        return []

    def set_mz_filter(self, mz: float):
        """Filter the results to measured features at mz (± precursor tolerance)."""
        try:
            tolerance: float = self.matching.precursor_tolerance_da(mz)
        except tk.TclError:
            messagebox.showerror("Invalid input", "Check the precursor tolerances.")
            return
        self.filters.set_mz_range(mz - tolerance, mz + tolerance)
        self.apply_filters()

    def on_click_mz_from_compound(self):
        if self.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        MzFromCompoundDialog(self)

    def on_click_remove_filters(self):
        """Empty all filter fields and show all matches."""
        self.filters.clear()
        self._show_filtered_matches(self.matches_table)
        self.set_status(f"Showing all {len(self.matches_table)} matches.")

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
        if (match := self.results.selected_match()) is None:
            messagebox.showerror("Error", "Please select a match row first.")
            return
        f_id_meas, f_id_lib = int(match.feature_id_meas), int(match.feature_id_lib)
        match_result: dict = next(m for m in self.matches[f_id_meas] if m['feature_id'] == f_id_lib)
        added: bool = self.add_pending_annotation(AnnotationStore.library_key(f_id_meas, f_id_lib), dict(
            kind="library",
            name=match_result.get('name'),
            formula=match_result.get('formula'),
            adduct=match_result.get('adduct'),
            ms2_score=match_result.get('ms2_score'),
            reference=f_id_lib,
            library=self.lib_db,
        ))
        if added:
            self.set_status(f"Added library feature {f_id_lib} for feature {f_id_meas} to the annotations. "
                            f"Use \"Manage annotations\" to write them to the database.")
        else:
            self.set_status(f"Library feature {f_id_lib} is already annotated for feature {f_id_meas}.")

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
        self.results.refresh_annotation_column()
        self.results.set_unsaved_annotation_count(sum(len(keys) for keys in self.annotations.unsaved_changes()))
        if self.annotation_window is not None:
            self.annotation_window.refresh()
        self.scatter.mark_annotated_features()

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
        """Row of the results table selected: mark its feature in the scatter plot and show the match."""
        f_id_meas, f_id_lib = int(match.feature_id_meas), int(match.feature_id_lib)
        self.scatter.mark_selected_feature(f_id_meas)
        if by_program:  # selected by the program (scatter click, restored selection): keep the overview
            return
        # e.g. selected again after sorting or filtering
        if ("match", (f_id_meas,), f_id_lib) != self.overview.shown_key:
            self.clear_comparison()
            self.overview.request("match", [f_id_meas], f_id_lib)

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
            match = self.results.selected_match()
            self.scatter.mark_selected_feature(f_id)
            self.overview.request("match", [f_id], int(match.feature_id_lib))
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
