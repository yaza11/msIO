"""Table of the library matches (sortable, copyable) and the annotation buttons."""

import tkinter as tk
from tkinter import ttk
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from msIO.gui.annotation.annotations import AnnotationStore
from msIO.gui.annotation.common import CollapsibleSection
from msIO.gui.annotation.matching import empty_matches_table

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI


class ResultsPanel(CollapsibleSection):
    """
    Filtered matches in a table:
        click on a header       sort ascending -> descending -> default order (Ctrl+click: further sort key)
        select a row            show the match in the overview (AnnotationGUI.on_match_selected)
        shift+click on a row    add its measured feature to the comparison
        right click / Ctrl+C    copy cell, row or table
    The buttons below collect annotations (see AnnotationGUI).
    """

    # column name in the matches table: (heading, width, anchor), this defines the order
    # "annotation" is not in the matches table, it shows the state of the manual annotation of the match
    COLUMNS: dict[str, tuple[str, int, str]] = {
        "feature_id_meas": ("Meas. ID", 60, "center"),
        "mz_meas": ("m/z", 80, "e"),
        "RT": ("RT (min)", 60, "e"),
        "ccs_meas": ("CCS (Å²)", 65, "e"),
        "intensity": ("Σ intensity", 80, "e"),
        "i1_i0": ("I1 / I0", 60, "e"),
        "name": ("Compound name", 180, "w"),
        "formula": ("Formula", 110, "w"),
        "adduct": ("Adduct", 80, "w"),
        "ms2_score": ("MS2 score", 70, "e"),
        "n_hits_ms2": ("#MS2 peaks", 75, "e"),
        "dmz_mda": ("Δm/z (mDa)", 80, "e"),
        "dmz_ppm": ("Δm/z (ppm)", 80, "e"),
        "dccs_percent": ("ΔCCS (%)", 70, "e"),
        "xlogp": ("logP", 55, "e"),
        "annotation": ("Annotation", 85, "center"),
        "source_library": ("Source", 150, "w"),
        "library": ("Library", 150, "w"),
    }
    # columns only shown if a row has a value (e.g. logP, which most libraries don't have)
    OPTIONAL_COLUMNS: tuple[str, ...] = ("xlogp",)
    # columns of the matches table shown instead (formulas that don't fit the measured m/z are in brackets)
    DISPLAY_COLUMNS: dict[str, str] = {"formula": "formula_display"}
    # background colors of the rows, alternating between measured features
    ROW_COLORS: list[str] = ["#E5FFCC", "#FFCCCC", "#CCFFCC", "#CCFFFF", "#CCCCFF", "#FFE5CC"]

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, "Results", stretch=True)
        self.app: "AnnotationGUI" = app
        self.columns: list[str] = list(self.COLUMNS)
        self.headings: dict[str, str] = {column: heading for column, (heading, _, _) in self.COLUMNS.items()}
        # shown (filtered and sorted) matches, the rows of the table are in this order
        self.table: pd.DataFrame = empty_matches_table(self.columns)
        # (column, descending) per clicked header, the first one is the primary key; empty: default order
        self.sort_keys: list[tuple[str, bool]] = []
        # set while selection events caused by the program are pending, see _ignore_pending_select_events
        self.ignore_select: bool = False

        tk_table_frame = ttk.Frame(self.body)
        tk_table_frame.pack(fill="both", expand=True)
        self.tk_tree: ttk.Treeview = ttk.Treeview(tk_table_frame, columns=self.columns, show="headings",
                                                  selectmode="browse")
        for column, (heading, width, anchor) in self.COLUMNS.items():
            self.tk_tree.heading(column, text=heading, command=lambda c=column: self.on_click_heading(c))
            self.tk_tree.column(column, width=width, anchor=anchor,
                                stretch=column in ("name", "source_library", "library"))
        # dummy rows only hold the properties of the measured feature
        self.tk_tree.tag_configure("dummy", foreground="gray40")
        self.tk_tree.bind("<<TreeviewSelect>>", self.on_select)
        self.tk_tree.bind("<ButtonRelease-1>", self.on_release)
        self.tk_tree.bind("<Control-ButtonRelease-1>", self.on_ctrl_release)
        self.tk_tree.bind("<Shift-Button-1>", self.on_shift_click)
        # the release of a shift+click must not show the match of the selected row (on_release)
        self.tk_tree.bind("<Shift-ButtonRelease-1>", lambda event: "break")
        # treeview cells can't be selected, so copying works with Ctrl+C (selected row) and a context menu
        self.tk_menu = tk.Menu(self.tk_tree, tearoff=False)
        self.tk_tree.bind("<Button-3>", self.on_right_click)
        self.tk_tree.bind("<Control-c>", lambda event: self.copy_selected_row())

        y_scroll = ttk.Scrollbar(tk_table_frame, orient="vertical", command=self.tk_tree.yview)
        x_scroll = ttk.Scrollbar(tk_table_frame, orient="horizontal", command=self.tk_tree.xview)
        self.tk_tree.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.tk_tree.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")
        tk_table_frame.rowconfigure(0, weight=1)
        tk_table_frame.columnconfigure(0, weight=1)

        tk_buttons = ttk.Frame(self.body)
        tk_buttons.pack(fill="x", pady=(5, 0))
        ttk.Button(tk_buttons, text="Export annotations as library ...",
                   command=app.on_click_export_library).pack(side="left")
        self.tk_manage_button = ttk.Button(tk_buttons, text="Manage annotations",
                                           command=app.on_click_manage_annotations)
        self.tk_manage_button.pack(side="right")
        ttk.Button(tk_buttons, text="Add library match", command=app.on_click_add_library_match).pack(
            side="right", padx=5)
        ttk.Button(tk_buttons, text="Add manual annotation", command=app.on_click_add_manual_annotation).pack(
            side="right")
        ttk.Button(tk_buttons, text="Mark fragment / duplicate / isotope / adduct ...",
                   command=app.on_click_mark_relation).pack(
            side="right", padx=5)

    # ---- content ----------------------------------------------------------------------------------------------------
    def set_table(self, table: pd.DataFrame):
        """Sort (by the clicked columns) and show matches, the selected match stays selected."""
        selected = self.selected_match()
        if not table.empty:
            table = self._sort(table)
        self.table = table.reset_index(drop=True)
        self._fill(self.table)

        if selected is not None and not self.table.empty:
            is_selected = ((self.table.feature_id_meas == selected.feature_id_meas)
                           & (self.table.feature_id_lib == selected.feature_id_lib)
                           & (self.table.library == selected.library))
            if is_selected.any():
                # restoring the selection must not change the overview (e.g. a feature shown from the scatter plot)
                self._ignore_pending_select_events()
                item: str = self.tk_tree.get_children()[int(np.argmax(is_selected))]
                self.tk_tree.selection_set(item)
                self.tk_tree.see(item)

    def clear(self):
        self.table = empty_matches_table(self.columns)
        self.tk_tree.delete(*self.tk_tree.get_children())

    def _fill(self, table: pd.DataFrame):
        """Replace the rows of the table, rows of the same measured feature share a background color."""
        self.tk_tree.delete(*self.tk_tree.get_children())
        # decided by all matches, so the columns don't come and go while filtering
        matches: pd.DataFrame = self.app.matches_table
        hidden: set[str] = {column for column in self.OPTIONAL_COLUMNS
                            if column not in matches.columns or matches[column].isna().all()}
        self.tk_tree.configure(displaycolumns=[column for column in self.columns if column not in hidden])
        if table.empty:
            return
        f_ids_meas = sorted(table['feature_id_meas'].unique())
        for color_idx, f_id_meas in enumerate(f_ids_meas):
            self.tk_tree.tag_configure(f"feature_{f_id_meas}",
                                       background=self.ROW_COLORS[color_idx % len(self.ROW_COLORS)])

        # values of all rows at once (iterrows takes seconds for tens of thousands of matches)
        source_columns: list[str] = [self.DISPLAY_COLUMNS.get(column, column) for column in self.columns]
        cells = pd.DataFrame({
            column: table[column] if column in table.columns else pd.Series("", index=table.index)
            for column in source_columns
        }).astype(object)
        rows: list[list] = cells.where(cells.notna(), "").to_numpy().tolist()
        idx_annotation: int = self.columns.index("annotation")
        annotations: AnnotationStore = self.app.annotations
        for values, row in zip(rows, table[['feature_id_meas', 'library', 'feature_id_lib', 'is_dummy']].itertuples()):
            values[idx_annotation] = self._annotation_state(row)
            tags = (f"feature_{row.feature_id_meas}", "dummy") if row.is_dummy else (f"feature_{row.feature_id_meas}",)
            self.tk_tree.insert("", 'end', values=values, tags=tags)

    def _annotation_state(self, row) -> str:
        """State of the annotation with the library match of a row (dummy rows have none)."""
        if row.is_dummy:
            return ""
        annotations: AnnotationStore = self.app.annotations
        return annotations.state(annotations.library_key(row.feature_id_meas, row.library, row.feature_id_lib))

    def refresh_annotation_column(self):
        for item, row in zip(self.tk_tree.get_children(), self.table.itertuples()):
            self.tk_tree.set(item, "annotation", self._annotation_state(row))

    def set_unsaved_annotation_count(self, n_changes: int):
        self.tk_manage_button.config(
            text=f"Manage annotations ({n_changes} unsaved)" if n_changes else "Manage annotations")

    # ---- sorting ----------------------------------------------------------------------------------------------------
    def _sort(self, table: pd.DataFrame) -> pd.DataFrame:
        """
        Sort matches by the sort keys, the first one is the primary key. Text is sorted case-insensitive, missing values
        last, ties keep the default order (as found by the search).
        """
        table = table.sort_values("default_order", kind="mergesort")
        # stable sorts from the last key to the first one give the sorting by several keys
        for column, descending in reversed(self.sort_keys):
            if column == "annotation":  # not a column of the matches table
                values = pd.Series([self._annotation_state(row) for row in table.itertuples()], index=table.index)
            else:
                values = table[column]
            if values.dtype == object:
                values = values.astype(str).str.lower()
            table = table.loc[values.sort_values(ascending=not descending, na_position="last", kind="mergesort").index]
        return table

    def on_click_heading(self, column: str, secondary: bool = False):
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
        for c, heading in self.headings.items():
            self.tk_tree.heading(c, text=heading + arrows.get(c, ""))
        self.set_table(self.table)

    def on_ctrl_release(self, event):
        """Ctrl+click on a column header: further sort key."""
        if self.tk_tree.identify_region(event.x, event.y) != "heading":
            return None
        column: str = self._column_at(event.x)
        self.on_click_heading(column, secondary=True)
        return "break"  # don't invoke the normal sorting of the header

    def _column_at(self, x: int) -> str:
        """Column at the x position (the shown columns can differ from self.columns, see OPTIONAL_COLUMNS)."""
        return self.tk_tree.column(self.tk_tree.identify_column(x), "id")

    # ---- selection --------------------------------------------------------------------------------------------------
    def selected_match(self) -> pd.Series | None:
        """Selected row of the shown matches."""
        selection: tuple[str, ...] = self.tk_tree.selection()
        if not selection:
            return None
        try:
            return self.table.iloc[self.tk_tree.index(selection[0])]
        except IndexError:
            return None

    def select_best_match(self, f_id_meas: int, show_match: bool = True) -> bool:
        """
        Select the match with the highest MS2 score of a measured feature (selecting a row shows the match in the
        overview, unless show_match is False). Returns False if the shown matches have none for the feature.
        """
        rows = self.table.index[self.table.feature_id_meas == f_id_meas]
        if rows.empty:
            return False
        # matches without MS2 score come last
        best_row = self.table.loc[rows, 'ms2_score'].fillna(-np.inf).idxmax()
        # rows of the table are in the order of self.table (its index is reset)
        item: str = self.tk_tree.get_children()[best_row]
        if not show_match:
            self._ignore_pending_select_events()
        self.tk_tree.selection_set(item)
        self.tk_tree.see(item)
        return True

    def _ignore_pending_select_events(self):
        """
        Selection events caused by the program (e.g. deleting and restoring the selected row) arrive later through the
        event queue, possibly several of them: ignore all of them until the queue is processed.
        """
        self.ignore_select = True
        self.app.tk_root.after_idle(lambda: setattr(self, "ignore_select", False))

    def on_select(self, event=None):
        # also fired with an empty selection when the table is refilled
        if (match := self.selected_match()) is not None:
            self.app.on_match_selected(match, by_program=self.ignore_select)

    def on_release(self, event):
        """Click on the already selected row (no selection event): show its match again, e.g. after a feature."""
        if (item := self.tk_tree.identify_row(event.y)) and item in self.tk_tree.selection():
            self.on_select()

    def on_shift_click(self, event):
        """Shift+click on a row: add its measured feature to the comparison (like shift+click in the scatter plot)."""
        if not (item := self.tk_tree.identify_row(event.y)):
            return None
        row = self.table.iloc[self.tk_tree.index(item)]
        self.app.add_feature_to_comparison(int(row.feature_id_meas))
        return "break"  # keep the selection

    # ---- copying ----------------------------------------------------------------------------------------------------
    def _copy_to_clipboard(self, text: str, description: str):
        self.app.tk_root.clipboard_clear()
        self.app.tk_root.clipboard_append(text)
        self.app.set_status(f"Copied {description} to the clipboard.")

    def _get_row_text(self, item: str) -> str:
        return "\t".join(str(value) for value in self.tk_tree.item(item, "values"))

    def copy_selected_row(self):
        if selection := self.tk_tree.selection():
            self._copy_to_clipboard(self._get_row_text(selection[0]), "the selected row")

    def copy_table(self):
        """Copy the shown (filtered and sorted) table as tab separated text (e.g. to paste into Excel)."""
        lines = ["\t".join(self.headings[column] for column in self.columns)]
        lines += [self._get_row_text(item) for item in self.tk_tree.get_children()]
        self._copy_to_clipboard("\n".join(lines), f"{len(lines) - 1} rows")

    def on_right_click(self, event):
        """Select the row under the mouse and offer to copy the cell, the row or the whole table."""
        item: str = self.tk_tree.identify_row(event.y)
        if not item:
            return
        self.tk_tree.selection_set(item)
        column: str = self._column_at(event.x)
        value: str = str(self.tk_tree.set(item, column))
        self.tk_menu.delete(0, "end")
        self.tk_menu.add_command(
            label=f"Copy {self.headings[column]}: {value if len(value) <= 40 else value[:37] + '...'}",
            command=lambda: self._copy_to_clipboard(value, f'"{value}"')
        )
        self.tk_menu.add_command(label="Copy row", command=self.copy_selected_row)
        self.tk_menu.add_command(label="Copy table", command=self.copy_table)
        self.tk_menu.tk_popup(event.x_root, event.y_root)
