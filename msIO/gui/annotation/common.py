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


def focused_widget(widget: tk.Misc) -> tk.Misc | None:
    """Widget with the keyboard focus in the application of widget, None if there is none (or tkinter can't resolve
    it, e.g. the open list of a combobox)."""
    try:
        return widget.focus_get()
    except KeyError:
        return None


def canvas_has_focus(canvas: FigureCanvasTkAgg | None) -> bool:
    return canvas is not None and focused_widget(canvas.get_tk_widget()) is canvas.get_tk_widget()


def create_canvas(fig, master: tk.Misc, focus_canvas: bool = False) -> FigureCanvasTkAgg:
    """
    Canvas of a figure that leaves the keyboard focus where it is. matplotlib focuses every new canvas, so redrawing
    the overview (on every selected row) would take the focus from the table: arrow keys and Ctrl+C would stop
    working there, as would the keys of the scatter plot. focus_canvas: give the new canvas the focus (e.g. it
    replaces a canvas that had it).
    """
    # setting the focus back afterwards is not enough: on Windows, Tk applies focus changes later, so the canvas would
    # still get it. matplotlib calls focus_set of the tk canvas in FigureCanvasTk.__init__, skip that call.
    original_focus_set = tk.Canvas.focus_set
    tk.Canvas.focus_set = lambda self: None
    try:
        canvas = FigureCanvasTkAgg(fig, master=master)
    finally:
        tk.Canvas.focus_set = original_focus_set
    if focus_canvas:
        canvas.get_tk_widget().focus_set()
    _delay_resize(canvas)
    return canvas


# a resized plot is redrawn once the size rests for this long (ms)
RESIZE_DELAY_MS: int = 150


def _delay_resize(canvas: FigureCanvasTkAgg):
    """
    matplotlib redraws the whole figure on every size change, which makes dragging a divider lag (one redraw per
    pixel). Redraw once the size rests instead, until then the old image is shown.
    """
    tk_widget = canvas.get_tk_widget()
    pending: dict = {"job": None, "event": None}

    def resize_later(event):
        pending["event"] = event
        if pending["job"] is not None:
            tk_widget.after_cancel(pending["job"])
        pending["job"] = tk_widget.after(RESIZE_DELAY_MS, resize_now)

    def resize_now():
        pending["job"] = None
        try:
            canvas.resize(pending["event"])
        except tk.TclError:  # destroyed meanwhile
            pass

    # replaces the binding of matplotlib (FigureCanvasTk binds <Configure> to resize)
    tk_widget.bind("<Configure>", resize_later)


class DividedArea(tk.PanedWindow):
    """
    Panes with movable dividers. While a divider is dragged only a line follows the mouse, the panes are resized once
    it is released: resizing them live re-lays out the whole window (tables, plots) on every mouse move, which makes
    dragging lag.
    """

    def __init__(self, master: tk.Misc, orient: str):
        background: str = ttk.Style(master).lookup("TFrame", "background") or "SystemButtonFace"
        super().__init__(master, orient=orient, opaqueresize=False, sashwidth=6, sashrelief="flat",
                         borderwidth=0, background=background)
        self.horizontal: bool = orient == "horizontal"

    def add_pane(self, widget: tk.Misc, grows: bool = True):
        """Add a pane, grows: it gets a share of additional space when the area grows."""
        self.add(widget, stretch="always" if grows else "never")

    def set_grows(self, widget: tk.Misc, grows: bool):
        self.paneconfigure(widget, stretch="always" if grows else "never")

    def divider_position(self, index: int = 0) -> int:
        """Position of a divider in pixels (x for horizontal, y for vertical areas)."""
        x, y = self.sash_coord(index)
        return x if self.horizontal else y

    def set_divider_position(self, position: int, index: int = 0):
        if self.horizontal:
            self.sash_place(index, position, 0)
        else:
            self.sash_place(index, 0, position)


class Tooltip:
    """Text shown in a small window while the mouse rests on a widget (after DELAY_MS)."""

    DELAY_MS: int = 600
    WRAP_LENGTH: int = 420

    def __init__(self, widget: tk.Misc, text: str):
        self.widget: tk.Misc = widget
        self.text: str = text
        self.job: str = None
        self.window: tk.Toplevel = None
        widget.bind("<Enter>", self._schedule, add="+")
        for sequence in ("<Leave>", "<ButtonPress>", "<KeyPress>"):
            widget.bind(sequence, self._hide, add="+")

    def _schedule(self, event=None):
        self._hide()
        self.job = self.widget.after(self.DELAY_MS, self._show)

    def _show(self):
        self.job = None
        if not self.text or not self.widget.winfo_viewable():
            return
        self.window = tk.Toplevel(self.widget)
        self.window.wm_overrideredirect(True)
        self.window.attributes("-topmost", True)
        tk.Label(self.window, text=self.text, justify="left", wraplength=self.WRAP_LENGTH, background="#FFFFE0",
                 relief="solid", borderwidth=1, padx=6, pady=4).pack()
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self.window.wm_geometry(f"+{x}+{y}")

    def _hide(self, event=None):
        if self.job is not None:
            self.widget.after_cancel(self.job)
            self.job = None
        if self.window is not None:
            self.window.destroy()
            self.window = None


def add_tooltip(text: str, *widgets: tk.Misc) -> None:
    """Show text while the mouse rests on any of the widgets (e.g. a label and its entry)."""
    for widget in widgets:
        Tooltip(widget, text)


def label_with_tooltip(master: tk.Misc, text: str, tooltip: str, *related_widgets: tk.Misc) -> ttk.Label:
    """Label whose tooltip is also shown for the related widgets (e.g. its entry)."""
    label = ttk.Label(master, text=text)
    add_tooltip(tooltip, label, *related_widgets)
    return label


def checkbutton_with_tooltip(master: tk.Misc, text: str, variable: tk.BooleanVar, command: Callable[[], None],
                             tooltip: str) -> ttk.Checkbutton:
    checkbutton = ttk.Checkbutton(master, text=text, variable=variable, command=command)
    add_tooltip(tooltip, checkbutton)
    return checkbutton


SortKeys = list[tuple[str, bool]]  # (column, descending), the first one is the primary key


def next_sort_keys(sort_keys: SortKeys, column: str, further_key: bool = False) -> SortKeys:
    """
    Sort keys after clicking a column header: ascending -> descending -> unsorted (default order). With further_key
    (Ctrl+click), the column is added to the keys (or cycled through the same states) instead of replacing them.
    """
    position = next((i for i, (sorted_column, _) in enumerate(sort_keys) if sorted_column == column), None)
    if not further_key or not sort_keys:
        if position == 0:  # primary key: reverse, then back to the default order
            return [(column, True)] + sort_keys[1:] if not sort_keys[0][1] else []
        return [(column, False)]
    if position is None:
        return sort_keys + [(column, False)]
    if not sort_keys[position][1]:
        return sort_keys[:position] + [(column, True)] + sort_keys[position + 1:]
    return sort_keys[:position] + sort_keys[position + 1:]


def show_sort_keys(tree: ttk.Treeview, headings: dict[str, str], sort_keys: SortKeys):
    """Arrows in the sorted column headers (numbered for several keys)."""
    numbered: bool = len(sort_keys) > 1
    arrows = {column: (" ▼" if descending else " ▲") + (str(i + 1) if numbered else "")
              for i, (column, descending) in enumerate(sort_keys)}
    for column, heading in headings.items():
        tree.heading(column, text=heading + arrows.get(column, ""))


class TreeSorter:
    """
    Sorting of a ttk.Treeview by clicking its column headers, as in the results table: ascending -> descending ->
    order of insertion, Ctrl+click adds a further sort key. Values are compared as numbers where possible (empty
    values last), otherwise as text (case-insensitive). Call apply() after refilling the table.
    """

    def __init__(self, tree: ttk.Treeview, headings: dict[str, str]):
        self.tree: ttk.Treeview = tree
        self.headings: dict[str, str] = headings
        self.keys: SortKeys = []
        for column in headings:
            tree.heading(column, command=lambda clicked=column: self.on_click(clicked))
        tree.bind("<Control-ButtonRelease-1>", self._on_ctrl_release, add="+")

    def _on_ctrl_release(self, event):
        if self.tree.identify_region(event.x, event.y) != "heading":
            return None
        self.on_click(self.tree.column(self.tree.identify_column(event.x), "id"), further_key=True)
        return "break"

    def on_click(self, column: str, further_key: bool = False):
        self.keys = next_sort_keys(self.keys, column, further_key)
        show_sort_keys(self.tree, self.headings, self.keys)
        self.apply()

    @staticmethod
    def _sort_value(text: str) -> tuple:
        """Numbers before text, empty values last (in both directions, see apply)."""
        try:
            return 0, float(text), ""
        except ValueError:
            return 1, 0., text.lower()

    def apply(self):
        """Order the rows by the sort keys (stable, ties keep the order of insertion)."""
        # order of insertion: the item ids (numbers given by the table, else Tk's I001, I002, ...)
        items: list[str] = sorted(self.tree.get_children(), key=self._insertion_position)
        # stable sorts from the last key to the first one give the sorting by several keys
        for column, descending in reversed(self.keys):
            filled = [item for item in items if str(self.tree.set(item, column)).strip()]
            empty = [item for item in items if not str(self.tree.set(item, column)).strip()]
            filled.sort(key=lambda item: self._sort_value(str(self.tree.set(item, column))), reverse=descending)
            items = filled + empty
        for position, item in enumerate(items):
            self.tree.move(item, "", position)

    @staticmethod
    def _insertion_position(item: str) -> tuple:
        return (0, int(item), "") if item.isdigit() else (1, 0, item)


def pack_canvas(canvas: FigureCanvasTkAgg):
    """Fill the parent with the canvas, the figure is scaled to the available space."""
    tk_widget = canvas.get_tk_widget()
    # without this, the canvas requests the full figure size and pushes the dividers of the panes
    tk_widget.config(width=1, height=1)
    tk_widget.pack(fill=tk.BOTH, expand=True)


class CollapsibleSection(ttk.Frame):
    """
    Section with a header line that collapses / expands its body (click the header). Widgets go into body.

    The body sits in a scrollable area: it asks for the size of its content and fills the space it gets, if it gets
    less (small window or pane), scroll bars appear (mouse wheel: vertical, Shift+wheel: horizontal; not over
    tables, lists and plots, which use the wheel themselves).

    stretch: the section fills the remaining space of its (packed) parent while it is expanded.
    on_toggle(expanded) is called after collapsing / expanding, e.g. to move the divider of a paned window.
    """

    # the size of the content is checked this often (ms), Tk has no event for a changed requested size
    SIZE_CHECK_MS: int = 250

    def __init__(self, master: tk.Misc, title: str, expanded: bool = True, stretch: bool = False,
                 on_toggle: Callable[[bool], None] = None, tooltip: str = None, **kwargs):
        super().__init__(master, relief="groove", borderwidth=1, padding=(2, 1), **kwargs)
        self.title: str = title
        self.expanded: bool = expanded
        self.stretch: bool = stretch
        self.on_toggle: Callable[[bool], None] | None = on_toggle

        self.tk_header = ttk.Label(self, cursor="hand2", font=("TkDefaultFont", 9, "bold"), padding=(4, 2))
        self.tk_header.pack(fill="x")
        self.tk_header.bind("<Button-1>", lambda event: self.toggle())
        add_tooltip((tooltip + "\n\n" if tooltip else "") + "Click the header to collapse / expand the section.",
                    self.tk_header)

        # body in a canvas, which scrolls it if it does not fit
        self.tk_viewport = ttk.Frame(self)
        self.tk_viewport.rowconfigure(0, weight=1)
        self.tk_viewport.columnconfigure(0, weight=1)
        background: str = ttk.Style(self).lookup("TFrame", "background") or "SystemButtonFace"
        self.tk_canvas = tk.Canvas(self.tk_viewport, highlightthickness=0, borderwidth=0, background=background)
        self.tk_canvas.grid(row=0, column=0, sticky="nsew")
        self.tk_y_scroll = ttk.Scrollbar(self.tk_viewport, orient="vertical", command=self.tk_canvas.yview)
        self.tk_x_scroll = ttk.Scrollbar(self.tk_viewport, orient="horizontal", command=self.tk_canvas.xview)
        self.tk_canvas.configure(yscrollcommand=self.tk_y_scroll.set, xscrollcommand=self.tk_x_scroll.set)
        self.body = ttk.Frame(self.tk_canvas, padding=(6, 2, 6, 4))
        self._body_window = self.tk_canvas.create_window(0, 0, window=self.body, anchor="nw")
        self._content_size: tuple[int, int] = None
        self.tk_canvas.bind("<Configure>", lambda event: self._fit_body())
        _bind_section_scrolling(self)
        self.after(self.SIZE_CHECK_MS, self._check_content_size)

        if expanded:
            self.tk_viewport.pack(fill="both", expand=True)
        self._update_header()

    # ---- scrolling --------------------------------------------------------------------------------------------------
    def _check_content_size(self):
        """The content changed its size (e.g. widgets added): ask for the new size and fit the body again."""
        try:
            size = self.body.winfo_reqwidth(), self.body.winfo_reqheight()
            if size != self._content_size:
                self._content_size = size
                # the section asks for the size of its content, as it would without scrolling
                self.tk_canvas.configure(width=size[0], height=size[1])
                self._fit_body()
            self.after(self.SIZE_CHECK_MS, self._check_content_size)
        except tk.TclError:  # destroyed
            pass

    def _fit_body(self):
        """Body as large as the visible area, but at least as large as its content; scroll bars for what is cut off."""
        if self._content_size is None:
            return
        content_width, content_height = self._content_size
        visible_width, visible_height = self.tk_canvas.winfo_width(), self.tk_canvas.winfo_height()
        if visible_width <= 1:  # not shown yet
            return
        width, height = max(visible_width, content_width), max(visible_height, content_height)
        self.tk_canvas.itemconfigure(self._body_window, width=width, height=height)
        self.tk_canvas.configure(scrollregion=(0, 0, width, height))
        # showing a scroll bar makes the visible area smaller, which calls this again
        for scrollbar, cut_off, sticky, position in (
                (self.tk_y_scroll, content_height > visible_height, "ns", dict(row=0, column=1)),
                (self.tk_x_scroll, content_width > visible_width, "ew", dict(row=1, column=0))):
            if cut_off and not scrollbar.winfo_ismapped():
                scrollbar.grid(sticky=sticky, **position)
            elif not cut_off and scrollbar.winfo_ismapped():
                scrollbar.grid_remove()

    def scroll(self, units: int, horizontal: bool = False) -> bool:
        """Scroll by units if the content is cut off in that direction, False if there is nothing to scroll."""
        scrollbar = self.tk_x_scroll if horizontal else self.tk_y_scroll
        if not scrollbar.winfo_ismapped():
            return False
        (self.tk_canvas.xview_scroll if horizontal else self.tk_canvas.yview_scroll)(units, "units")
        return True

    # ---- header -----------------------------------------------------------------------------------------------------
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
            self.tk_viewport.pack(fill="both", expand=True)
        else:
            self.tk_viewport.pack_forget()
        if self.stretch and self.winfo_manager() == "pack":
            self.pack_configure(expand=expanded)
        self._update_header()
        if self.on_toggle is not None:
            self.on_toggle(expanded)


# widgets that use the mouse wheel themselves (a tk.Canvas that is not the canvas of a section: matplotlib plot)
_OWN_WHEEL_SCROLLING = (ttk.Treeview, tk.Listbox, tk.Text, ttk.Combobox, ttk.Spinbox)


def _bind_section_scrolling(section: CollapsibleSection):
    """Mouse wheel over a section scrolls it (if its content is cut off), bound once for the application."""
    root: tk.Misc = section._root()
    if getattr(root, "_section_wheel_bound", False):
        return
    root._section_wheel_bound = True

    def on_wheel(event, horizontal: bool):
        widget = event.widget
        if not isinstance(widget, tk.Misc) or isinstance(widget, _OWN_WHEEL_SCROLLING):
            return None
        # the section around the widget under the mouse
        while widget is not None and not isinstance(widget, CollapsibleSection):
            if isinstance(widget, tk.Canvas) and not isinstance(widget.master.master, CollapsibleSection):
                return None  # a plot
            widget = widget.master
        if widget is not None:
            widget.scroll(-1 if event.delta > 0 else 1, horizontal)
        return None

    root.bind_all("<MouseWheel>", lambda event: on_wheel(event, False), add="+")
    root.bind_all("<Shift-MouseWheel>", lambda event: on_wheel(event, True), add="+")
