"""Constants, parsing of user input and small widgets shared by the parts of the annotation GUI."""

import re
import tkinter as tk
from tkinter import messagebox, ttk
from typing import Callable

import numpy as np
from LipidCalculator import Adduct, CompoundDict
from LipidCalculator.compound_creation.formula_parser import parse_equation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

try:
    from LipidCalculator.gui.common import ADDUCTS
except ImportError:  # LipidCalculator without GUI
    ADDUCTS: list[str] = ['[M+H]+', '[M+NH4]+', '[M+Na]+', '[M+K]+', '[M-H]-', '[M+Cl]-', '[M+HCOO]-']

SQLITE_FILE_TYPES: list[tuple[str, str]] = [("SQLite Files", "*.sqlite *.db *.sql"), ("All files", "*.*")]

# (kind, measured feature id, reference) identifies a manual annotation, kind is one of the
# ManualAnnotationKind of FeatureManagerDB, the reference is the library feature id ("library"),
# the reference feature id (relations to another feature, see FEATURE_RELATIONS) or "name|adduct" ("manual")
AnnotationKey = tuple[str, int, str]
# relations between measured features: kind -> (label, role of the feature, role of the reference feature)
FEATURE_RELATIONS: dict[str, tuple[str, str, str]] = {
    "in_source_fragment": ("in-source fragment", "Fragment", "Parent"),
    "duplicate": ("duplicate", "Duplicate", "Original"),
    "isotope": ("isotope peak", "Isotope peak", "Monoisotopic peak"),
    "adduct": ("adduct", "Adduct", "Main adduct"),
}
FEATURE_RELATION_KINDS: tuple[str, ...] = tuple(FEATURE_RELATIONS)
ANNOTATION_KIND_LABELS: dict[str, str] = {
    "library": "library match", "manual": "manual", **{kind: label for kind, (label, _, _) in FEATURE_RELATIONS.items()}
}

# colors of the two features compared by shift-click, used in the scatter plot and the overview
COMPARISON_COLORS: tuple[str, str] = ("C2", "C1")


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


# inputs of the IPL builder when it was last closed or used (session only), shared by all dialogs that open it
ipl_builder_state: dict = {}
# variables of SingleIPLFrame in the order they are restored, the fields of the core are reset when it changes, so
# they come after it (_update_chain_state); older LipidCalculator versions don't have all of them
IPL_BUILDER_VARIABLES: tuple[str, ...] = ("core_class_var", "head_var", "head2_var", "core_var", "chains_var",
                                          "split_var", "sort_var")
IPL_BUILDER_CORE_VARIABLES: tuple[str, ...] = ("core_chains_var", "core_db_var", "core_oh_var", "core_me_var")


def save_ipl_builder_state(frame):
    """Remember the inputs and the selected adduct of the IPL builder."""
    for name in (*IPL_BUILDER_VARIABLES, "rings_var", "isomer_var", *IPL_BUILDER_CORE_VARIABLES):
        if hasattr(frame, name):
            try:
                ipl_builder_state[name] = getattr(frame, name).get()
            except tk.TclError:  # e.g. empty number field
                pass
    row = frame.table.selected_row()
    ipl_builder_state["adduct"] = None if row is None else row["Adduct"]


def restore_ipl_builder_state(frame):
    """Fill the IPL builder with the remembered inputs and select the remembered adduct once it is built."""
    if not ipl_builder_state:
        return

    def restore(names):
        for name in names:
            if name in ipl_builder_state and hasattr(frame, name):
                getattr(frame, name).set(ipl_builder_state[name])

    restore(IPL_BUILDER_VARIABLES)
    if hasattr(frame, "_set_core_class"):  # cores of the remembered class in the list
        frame.core_cb.configure(values=_cores_of_class(frame.core_class_var.get()))
    frame._update_chain_state()  # shows the fields of the core and resets them
    if ipl_builder_state.get("rings_var") and hasattr(frame, "_update_isomers"):  # tetraethers
        restore(["rings_var"])
        frame._update_isomers()
        restore(["isomer_var"])
    restore(IPL_BUILDER_CORE_VARIABLES)

    def select_adduct(tries: int = 50):
        """The structure is built shortly after the inputs change, then its table has rows."""
        if not frame.winfo_exists():
            return
        tree = frame.table.tree
        if not tree.get_children():
            if tries:
                frame.after(100, lambda: select_adduct(tries - 1))
            return
        for item in tree.get_children():
            if tree.set(item, "Adduct") == ipl_builder_state.get("adduct"):
                tree.selection_set(item)
                tree.see(item)

    frame.after(100, select_adduct)


def _cores_of_class(core_class: str) -> list[str]:
    """Cores of a class in the IPL builder."""
    from LipidCalculator.gui.ipl_gui import CORE_CLASSES
    return CORE_CLASSES[core_class]


def open_ipl_builder(parent: tk.Misc, send_label: str, on_send: Callable[[object, dict, tk.Toplevel], None]):
    """
    Open the IPL builder of LipidCalculator in a window. Its button send_label calls on_send(frame, payload, window),
    the payload holds the neutral formula and the adduct selected in the table of the builder (if any).
    The inputs from the last time the builder was used (in this session) are restored.
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

    def send(payload):
        save_ipl_builder_state(frames[0])
        on_send(frames[0], payload, tk_window)

    def close():
        save_ipl_builder_state(frames[0])
        tk_window.destroy()

    frame = SingleIPLFrame(tk_window, send_targets={send_label: send})
    frames.append(frame)
    frame.pack(fill="both", expand=True)
    restore_ipl_builder_state(frame)
    tk_window.protocol("WM_DELETE_WINDOW", close)


# the batch window is hidden when closed and shown again, so its inputs and generated table are kept (session only)
_ipl_batch_window: tk.Toplevel | None = None


def batch_targets(frame, selected_only: bool) -> list:
    """
    m/z targets of the table of an IPL batch (BatchIPLFrame of LipidCalculator): one per row (the selected ones or
    all) and adduct column.
    """
    from msIO.gui.annotation.matching import MzTarget
    table = frame.table
    columns: list[str] = list(table.columns)
    adducts: list[str] = [column for column in columns if column in ADDUCTS or column.startswith("[")]
    targets: list[MzTarget] = []
    for values in table.rows(selected_only=selected_only):
        row: dict = dict(zip(columns, values))
        formula = re.sub(r"[+-]\d*$", "", row.get("Formula") or "") or None  # rdkit formulas can carry a charge
        for adduct in adducts:
            try:
                mz = float(row[adduct])
            except (TypeError, ValueError):
                continue
            targets.append(MzTarget(name=row.get("Name") or "?", formula=formula, adduct=adduct, mz=mz,
                                    smiles=row.get("SMILES") or None))
    return targets


def open_ipl_batch(root: tk.Misc, on_use: Callable[[list], None]):
    """
    Open the batch generation of LipidCalculator (e.g. archaeols with ranges of chain lengths and double bonds). Its
    button "Use for m/z filter" calls on_use(targets) with the targets of the selected rows (all rows if none is
    selected), see batch_targets. The window keeps its inputs and table until the session ends.
    """
    global _ipl_batch_window
    try:
        from LipidCalculator.gui.ipl_gui import BatchIPLFrame
    except Exception as error:
        messagebox.showerror("IPL batch", f"The batch generation of LipidCalculator is not available:\n{error}",
                             parent=root)
        return
    if _ipl_batch_window is not None and _ipl_batch_window.winfo_exists():
        _ipl_batch_window.on_use = on_use
        _ipl_batch_window.deiconify()
        _ipl_batch_window.lift()
        return _ipl_batch_window

    tk_window = tk.Toplevel(root)
    tk_window.title("IPL batch (generate, select rows or none for all, then \"Use for m/z filter\")")
    tk_window.geometry("1250x780")
    tk_window.on_use = on_use
    frame = BatchIPLFrame(tk_window)
    tk_window.frame = frame
    tk_buttons = ttk.Frame(tk_window, padding=(8, 0, 8, 8))
    tk_buttons.pack(side="bottom", fill="x")
    frame.pack(fill="both", expand=True)

    def use():
        selected: bool = bool(frame.table.tree.selection())
        if not (targets := batch_targets(frame, selected_only=selected)):
            messagebox.showerror("IPL batch", "Generate IPLs with at least one adduct first.", parent=tk_window)
            return
        tk_window.on_use(targets)

    ttk.Button(tk_buttons, text="Close", command=tk_window.withdraw).pack(side="right")
    ttk.Button(tk_buttons, text="Use for m/z filter", command=use).pack(side="right", padx=5)
    ttk.Label(tk_buttons, foreground="gray",
              text="m/z of each generated IPL and adduct ± precursor tolerance; features without library matches "
                   "are shown as dummy rows").pack(side="left")
    tk_window.protocol("WM_DELETE_WINDOW", tk_window.withdraw)
    _ipl_batch_window = tk_window
    return tk_window


def pack_canvas(canvas: FigureCanvasTkAgg):
    """Fill the parent with the canvas, the figure is scaled to the available space."""
    tk_widget = canvas.get_tk_widget()
    # without this, the canvas requests the full figure size and pushes the dividers of the panes
    tk_widget.config(width=1, height=1)
    tk_widget.pack(fill=tk.BOTH, expand=True)


class CollapsibleSection(ttk.Frame):
    """
    Section with a header line that collapses / expands its body (click the header). Widgets go into body.

    stretch: the section fills the remaining space of its (packed) parent while it is expanded.
    on_toggle(expanded) is called after collapsing / expanding, e.g. to move the divider of a paned window.
    """

    def __init__(self, master: tk.Misc, title: str, expanded: bool = True, stretch: bool = False,
                 on_toggle: Callable[[bool], None] = None, **kwargs):
        super().__init__(master, relief="groove", borderwidth=1, padding=(2, 1), **kwargs)
        self.title: str = title
        self.expanded: bool = expanded
        self.stretch: bool = stretch
        self.on_toggle: Callable[[bool], None] | None = on_toggle

        self.tk_header = ttk.Label(self, cursor="hand2", font=("TkDefaultFont", 9, "bold"), padding=(4, 2))
        self.tk_header.pack(fill="x")
        self.tk_header.bind("<Button-1>", lambda event: self.toggle())
        self.body = ttk.Frame(self, padding=(6, 2, 6, 4))
        if expanded:
            self.body.pack(fill="both", expand=True)
        self._update_header()

    def _update_header(self):
        self.tk_header.config(text=f"{'▼' if self.expanded else '▶'}  {self.title}")

    def set_title(self, title: str):
        self.title = title
        self._update_header()

    def header_height(self) -> int:
        """Height of the collapsed section in pixels."""
        return self.tk_header.winfo_reqheight() + 6

    def toggle(self):
        self.set_expanded(not self.expanded)

    def set_expanded(self, expanded: bool):
        if expanded == self.expanded:
            return
        self.expanded = expanded
        if expanded:
            self.body.pack(fill="both", expand=True)
        else:
            self.body.pack_forget()
        if self.stretch and self.winfo_manager() == "pack":
            self.pack_configure(expand=expanded)
        self._update_header()
        if self.on_toggle is not None:
            self.on_toggle(expanded)
