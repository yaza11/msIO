"""Constants, parsing of user input and small widgets shared by the parts of the annotation GUI."""

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
# the parent feature id ("in_source_fragment") or "name|adduct" ("manual")
AnnotationKey = tuple[str, int, str]
ANNOTATION_KIND_LABELS: dict[str, str] = {
    "library": "library match", "manual": "manual", "in_source_fragment": "in-source fragment"
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
