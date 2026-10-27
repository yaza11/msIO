"""Windows and dialogs of the annotation GUI."""

import threading
import tkinter as tk
from tkinter import messagebox, ttk
from typing import TYPE_CHECKING

import pandas as pd
from LipidCalculator import Adduct, CompoundDict

from msIO.gui.annotation.common import ADDUCTS, ANNOTATION_KIND_LABELS, AnnotationKey, open_ipl_builder

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI


class AnnotationManagerWindow(tk.Toplevel):
    """
    Window to review manual annotations before writing them to the measured database.

    Lists new annotations (library matches, manual annotations, in-source fragments) and the ones
    already in the database. Each can be kept or discarded; "Write to database" adds the kept new
    ones and removes the discarded ones from the database. The state lives in AnnotationGUI.annotations,
    so closing this window doesn't lose any changes.
    """

    def __init__(self, app: "AnnotationGUI"):
        super().__init__(app.tk_root)
        self.app: "AnnotationGUI" = app
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
        self.tk_table.bind("<space>", lambda event: self.app.toggle_annotations(self._get_selected_keys()))

        y_scroll = ttk.Scrollbar(tk_table_frame, orient="vertical", command=self.tk_table.yview)
        self.tk_table.configure(yscrollcommand=y_scroll.set)
        self.tk_table.pack(side="left", fill="both", expand=True)
        y_scroll.pack(side="left", fill="y")

        tk_buttons = ttk.Frame(self, padding=10)
        tk_buttons.pack(fill="x")
        ttk.Button(tk_buttons, text="Keep selected",
                   command=lambda: self.app.toggle_annotations(self._get_selected_keys(), keep=True)).pack(side="left")
        ttk.Button(tk_buttons, text="Discard selected",
                   command=lambda: self.app.toggle_annotations(self._get_selected_keys(), keep=False)).pack(
            side="left", padx=5)
        ttk.Button(tk_buttons, text="Keep all",
                   command=lambda: self.app.toggle_annotations(list(self.iid_to_key.values()), keep=True)).pack(
            side="left")
        ttk.Button(tk_buttons, text="Discard all",
                   command=lambda: self.app.toggle_annotations(list(self.iid_to_key.values()), keep=False)).pack(
            side="left", padx=5)
        ttk.Button(tk_buttons, text="Close", command=self.on_close).pack(side="right")
        ttk.Button(tk_buttons, text="Write to database",
                   command=lambda: self.app.write_annotations_to_db(parent=self)).pack(side="right", padx=5)

        self.refresh()

    def _get_selected_keys(self) -> list[AnnotationKey]:
        return [self.iid_to_key[iid] for iid in self.tk_table.selection()]

    def refresh(self):
        """Rebuild the rows from the annotation state, keeping selection and scroll position."""
        selected_keys = set(self._get_selected_keys())
        y_view = self.tk_table.yview()[0]
        self.tk_table.delete(*self.tk_table.get_children())
        self.iid_to_key = {}
        annotations = self.app.annotations

        # new annotations first, then the ones in the database, each sorted by measured feature
        def sort_key(item) -> tuple:
            (kind, f_id_meas, ref), _ = item
            return f_id_meas, kind, ref

        rows = ([(key, info, False) for key, info in sorted(annotations.pending.items(), key=sort_key)]
                + [(key, info, True) for key, info in sorted(annotations.saved.items(), key=sort_key)])
        for idx, (key, info, is_saved) in enumerate(rows):
            kind, f_id_meas, _ = key
            discarded: bool = key in annotations.discarded
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
            self.app.toggle_annotations([self.iid_to_key[iid]])

    def on_close(self):
        self.app.annotation_window = None
        self.destroy()


class ManualAnnotationDialog(tk.Toplevel):
    """
    Dialog to annotate a measured feature with a compound entered by hand.

    Name, formula, adduct and SMILES can be typed in or taken from the IPL builder of LipidCalculator.
    The m/z calculated from formula and adduct is compared to the m/z of the feature.
    """

    def __init__(self, app: "AnnotationGUI", f_id_meas: int | None):
        super().__init__(app.tk_root)
        self.app: "AnnotationGUI" = app
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
            mz_meas = self.app.meas_db.mzs.get(int(self.tk_f_id_meas.get()))
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
        if f_id_meas not in self.app.meas_db.mzs:
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

        added: bool = self.app.add_pending_annotation(("manual", f_id_meas, f"{name}|{adduct}"), dict(
            kind="manual", name=name, formula=formula, adduct=adduct, smiles=smiles, reference=None
        ))
        self.app.set_status(f"Added manual annotation {name} {adduct or ''} for feature {f_id_meas}." if added
                            else f"Feature {f_id_meas} is already annotated as {name} {adduct or ''}.")
        self.destroy()


class MzFromCompoundDialog(tk.Toplevel):
    """
    Dialog to set the m/z filter of the results table from a compound: searched by name in the library, built with
    the IPL builder of LipidCalculator or entered as formula and adduct. The filter is m/z ± precursor tolerance.
    """

    MAX_SEARCH_RESULTS: int = 500

    def __init__(self, app: "AnnotationGUI"):
        super().__init__(app.tk_root)
        self.app: "AnnotationGUI" = app
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

        if self.app.lib_db is None:
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
            tolerance = self.app.matching.precursor_tolerance_da(mz)
        except (ValueError, tk.TclError):
            self.tk_range_info.set("")
            return
        self.tk_range_info.set(f"filter: {mz - tolerance:.4f} - {mz + tolerance:.4f} (± {tolerance * 1e3:.2f} mDa, "
                               f"precursor tolerance)")

    def on_click_search(self):
        """Search the library names in a worker thread (millions of names in large libraries)."""
        if not (name := self.tk_name.get().strip()) or self.app.lib_db is None:
            return
        self.tk_search_button.state(["disabled"])
        self.tk_search_info.set("searching ...")
        lib_db, exact = self.app.lib_db, self.tk_exact.get()

        def worker():
            try:
                f_ids = lib_db.find_by_name(name, substring=not exact)
                rows = [(f_id, lib_db.names.get(f_id), lib_db.formula_metaboscape.get(f_id), lib_db.adduct.get(f_id),
                         lib_db.mzs.get(f_id)) for f_id in f_ids[:self.MAX_SEARCH_RESULTS]]
                result = (len(f_ids), rows, None)
            except Exception as error:
                result = (0, [], f"{type(error).__name__}: {error}")
            self.app.call_in_main_thread(lambda: self._show_search_results(*result))

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
        self.app.set_mz_filter(mz)


class InSourceFragmentDialog(tk.Toplevel):
    """
    Dialog to annotate a measured feature as in-source fragment of another one.

    Shows the differences in m/z and RT of both features and the marked mass differences connecting them.
    """

    def __init__(self, app: "AnnotationGUI", f_id_fragment: int | None, f_id_parent: int | None):
        super().__init__(app.tk_root)
        self.app: "AnnotationGUI" = app
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
        return f_ids if all(f_id in self.app.meas_db.mzs for f_id in f_ids) else None

    def _update_info(self):
        """Show Δm/z, ΔRT and the marked mass differences connecting both features, suggest the note from them."""
        if (f_ids := self._read_feature_ids()) is None:
            self.tk_info.set("Enter two measured feature IDs (shift-click both in the scatter plot to prefill).")
            return
        (rt_fragment, mz_fragment), (rt_parent, mz_parent) = (self.app.get_rt_mz(f_id) for f_id in f_ids)
        lines = [f"Δm/z (parent - fragment) {mz_parent - mz_fragment:.4f}, ΔRT {rt_parent - rt_fragment:+.2f} min"]
        if mz_parent < mz_fragment:
            lines.append("The parent has a lower m/z than the fragment.")
        if marked := self.app.scatter.find_marked_mass_differences(*f_ids):
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
        added: bool = self.app.add_pending_annotation(("in_source_fragment", f_id_fragment, str(f_id_parent)), dict(
            kind="in_source_fragment",
            name=f"in-source fragment of {f_id_parent}" + (f" ({note})" if note else ""),
            note=note,
            reference=f_id_parent,
        ))
        self.app.set_status(f"Marked feature {f_id_fragment} as in-source fragment of {f_id_parent}." if added
                            else f"Feature {f_id_fragment} is already marked as in-source fragment of {f_id_parent}.")
        self.destroy()
