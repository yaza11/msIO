"""
Keyboard shortcuts of the annotation GUI and the help window (F1) listing them together with the mouse actions.

The shortcuts are bound to the window of the GUI (dialogs have their own keys). They only act while the GUI is shown,
e.g. not while another tab of msIO.gui.main is selected. Shortcuts with Ctrl also work while typing in a field (Tk's
text fields use Ctrl+D, H, K and T, so these are left out), plain keys only elsewhere. Keys over the scatter plot are
handled by the plot (scatter_panel.PLOT_KEYS) and only listed here.
"""
import tkinter as tk
from dataclasses import dataclass
from tkinter import ttk
from typing import TYPE_CHECKING, Callable

from msIO.gui.annotation.scatter_panel import PLOT_KEYS

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI


@dataclass(frozen=True)
class Shortcut:
    keys: str  # as shown in the help
    description: str
    sequences: tuple[str, ...] = ()  # tk event sequences, empty: handled elsewhere and only listed
    action: Callable[["AnnotationGUI"], None] = None
    while_typing: bool = True  # also act while the focus is in a text field


def _focus_filters(app: "AnnotationGUI"):
    app.filters.set_expanded(True)
    app.filters.tk_name.focus_set()
    app.filters.tk_name.select_range(0, tk.END)


def _stop_comparing(app: "AnnotationGUI"):
    if app.compared_f_ids:
        app.clear_comparison()
        app.set_status("Stopped comparing features.")


# (group, shortcuts) in the order of the help window
SHORTCUTS: list[tuple[str, list[Shortcut]]] = [
    ("General", [
        Shortcut("F1", "Show the shortcuts", ("<F1>",), lambda app: show_shortcuts(app)),
        Shortcut("Ctrl+O", "Open the measured database", ("<Control-o>",), lambda app: app.on_browse_meas_db()),
        Shortcut("Ctrl+Shift+O", "Add libraries", ("<Control-O>",), lambda app: app.on_add_lib_db()),
        Shortcut("Ctrl+Enter", "Find matches", ("<Control-Return>", "<Control-KP_Enter>"),
                 lambda app: app.on_click_find_matches()),
    ]),
    ("Results table and filters", [
        Shortcut("Ctrl+F", "Go to the filters (name field)", ("<Control-f>",), _focus_filters),
        Shortcut("Ctrl+Shift+F", "Remove all filters", ("<Control-F>",), lambda app: app.on_click_remove_filters()),
        Shortcut("Ctrl+T", "Go to the table (then ↑ / ↓: previous / next row)", ("<Control-t>",),
                 lambda app: app.results.focus_table()),
        Shortcut("Ctrl+↓ / Ctrl+↑", "Next / previous feature in the table (skips further matches of a feature)",
                 ("<Control-Down>",), lambda app: app.results.select_adjacent_feature(1)),
        Shortcut("", "", ("<Control-Up>",), lambda app: app.results.select_adjacent_feature(-1)),
        Shortcut("Ctrl+C", "Copy the selected row (in the table)"),
    ]),
    ("Overview", [
        Shortcut("Ctrl+G", "Go to feature: type one ID (select it) or two (compare them), then Enter",
                 ("<Control-g>",), lambda app: app.overview.focus_feature_ids()),
        Shortcut("Esc", "Stop comparing features", ("<Escape>",), _stop_comparing, while_typing=False),
    ]),
    ("Annotations", [
        Shortcut("Ctrl+L", "Add the selected library match", ("<Control-l>",),
                 lambda app: app.on_click_add_library_match()),
        Shortcut("Ctrl+N", "Add a manual annotation", ("<Control-n>",),
                 lambda app: app.on_click_add_manual_annotation()),
        Shortcut("Ctrl+R", "Mark fragment / duplicate / isotope / adduct", ("<Control-r>",),
                 lambda app: app.on_click_mark_relation()),
        Shortcut("Ctrl+M", "Manage annotations", ("<Control-m>",), lambda app: app.on_click_manage_annotations()),
        Shortcut("Ctrl+S", "Write the annotation changes to the measured database", ("<Control-s>",),
                 lambda app: app.write_annotations_to_db()),
        Shortcut("Ctrl+E", "Export annotations as library", ("<Control-e>",), lambda app: app.on_click_export_library()),
    ]),
    ("Scatter plot (mouse over the plot)", [Shortcut(keys, description) for _, keys, description, _ in PLOT_KEYS]),
]

# (where, mouse action, effect)
MOUSE_ACTIONS: list[tuple[str, str, str]] = [
    ("Scatter plot", "Click", "Select the feature: its best match in the table, otherwise the feature alone"),
    ("Scatter plot", "Shift+click", "Compare two features (the older one is replaced)"),
    ("Scatter plot", "Ctrl+click", "Add to / remove from a series (after releasing Ctrl: new series)"),
    ("Results table", "Click a row", "Show the match (the plot pans to the feature if needed)"),
    ("Results table", "Shift+click a row", "Add its feature to the comparison"),
    ("Results table", "Click a header", "Sort: ascending → descending → default order"),
    ("Results table", "Ctrl+click a header", "Further sort key"),
    ("Results table", "Right click", "Copy the cell, the row or the table"),
    ("Overview", "Click an MS2 peak", "Annotate the peak"),
    ("Everywhere", "Click a section header", "Collapse / expand the section"),
]

_TEXT_FIELDS = (tk.Entry, ttk.Entry, tk.Text, tk.Spinbox)


def bind_shortcuts(app: "AnnotationGUI"):
    """Bind the shortcuts to the window of the GUI."""
    window: tk.Misc = app.tk_root.winfo_toplevel()

    def run(event, shortcut: Shortcut):
        if not app.tk_root.winfo_ismapped():  # e.g. another tab of the main window is shown
            return None
        if not shortcut.while_typing and isinstance(event.widget, _TEXT_FIELDS):
            return None
        shortcut.action(app)
        return "break"

    for _, shortcuts in SHORTCUTS:
        for shortcut in shortcuts:
            for sequence in shortcut.sequences:
                window.bind(sequence, lambda event, s=shortcut: run(event, s), add="+")


def show_shortcuts(app: "AnnotationGUI"):
    """Window listing the shortcuts and mouse actions (only one is opened, F1 / Esc close it again)."""
    existing: tk.Toplevel | None = getattr(app, "_shortcuts_window", None)
    if existing is not None and existing.winfo_exists():
        existing.destroy()  # F1 toggles the help
        return

    window = tk.Toplevel(app.tk_root)
    app._shortcuts_window = window
    window.title("Shortcuts and mouse actions")
    window.geometry("860x680")
    for key in ("<Escape>", "<F1>"):
        window.bind(key, lambda event: window.destroy())

    tree = ttk.Treeview(window, columns=("keys", "description"), show="tree headings")
    tree.heading("#0", text="")
    tree.heading("keys", text="Keys / mouse", anchor="w")
    tree.heading("description", text="Action", anchor="w")
    tree.column("#0", width=20, stretch=False)
    tree.column("keys", width=260, stretch=False)
    tree.column("description", width=540)
    tree.tag_configure("group", font=("TkDefaultFont", 9, "bold"))
    scrollbar = ttk.Scrollbar(window, orient="vertical", command=tree.yview)
    tree.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side="right", fill="y")
    tree.pack(fill="both", expand=True, padx=(8, 0), pady=8)

    for group, shortcuts in SHORTCUTS:
        parent = tree.insert("", "end", values=(group, ""), tags=("group",), open=True)
        for shortcut in shortcuts:
            if shortcut.keys:  # bindings without own entry are described by the previous one
                tree.insert(parent, "end", values=(shortcut.keys, shortcut.description))
    groups: dict[str, str] = {}
    for where, action, effect in MOUSE_ACTIONS:
        if where not in groups:
            groups[where] = tree.insert("", "end", values=(f"Mouse: {where}", ""), tags=("group",), open=True)
        tree.insert(groups[where], "end", values=(action, effect))
    window.focus_set()
