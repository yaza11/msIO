"""
Keyboard shortcuts of the annotation GUI and the help window (F1) listing them together with the mouse actions.

The shortcuts are bound to the window of the GUI (dialogs have their own keys). They only act while the GUI is shown,
e.g. not while another tab of msIO.gui.main is selected. Shortcuts with Ctrl also work while typing in a field (Tk's
text fields use Ctrl+D, H, K and T, so these are left out), plain keys only elsewhere. Keys over the scatter plot are
handled by the plot (scatter_panel.PLOT_KEYS) and only listed here. Esc closes pop-up windows (all windows of the
application except the main window, see close_popups_with_escape).
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
    only_while_typing: bool = False  # only act while the focus is in a text field


def _focus_filters(app: "AnnotationGUI"):
    app.filters.set_expanded(True)
    app.filters.tk_name.focus_set()
    app.filters.tk_name.select_range(0, tk.END)


def _leave_text_field(app: "AnnotationGUI"):
    """Esc in a text field: give the focus to the window, so that plain keys (V, M, Z, Esc) work again."""
    app.tk_root.winfo_toplevel().focus_set()


def _stop_comparing(app: "AnnotationGUI"):
    if app.compared_f_ids:
        app.clear_comparison()
        app.set_status("Stopped comparing features.")


# (group, shortcuts) in the order of the help window
SHORTCUTS: list[tuple[str, list[Shortcut]]] = [
    ("General", [
        Shortcut("F1 / ?", "Show the shortcuts (? button at the top right)", ("<F1>",), lambda app: show_shortcuts(app)),
        Shortcut("Esc", "In a text field: leave it (then plain keys like V, M, Z work again)", ("<Escape>",),
                 _leave_text_field, only_while_typing=True),
        Shortcut("Esc", "Close a pop-up window"),
        Shortcut("Ctrl+O", "Open the measured database", ("<Control-o>",), lambda app: app.on_browse_meas_db()),
        Shortcut("Ctrl+Shift+O", "Add libraries", ("<Control-O>",), lambda app: app.on_add_lib_db()),
        Shortcut("Ctrl+Enter", "Find matches", ("<Control-Return>", "<Control-KP_Enter>"),
                 lambda app: app.on_click_find_matches()),
    ]),
    ("Results table and filters", [
        Shortcut("Ctrl+F", "Go to the filters (name field)", ("<Control-f>",), _focus_filters),
        Shortcut("Ctrl+Shift+F", "Remove all filters", ("<Control-F>",), lambda app: app.on_click_remove_filters()),
        Shortcut("Ctrl+T", "Go to the table", ("<Control-t>",), lambda app: app.results.focus_table()),
        Shortcut("↑ / ↓", "Previous / next row in the table (after clicking it or Ctrl+T)"),
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
    ("Scatter plot", [
        Shortcut("Z", "Zoom to the features in the table", ("<Key-z>", "<Key-Z>"),
                 lambda app: app.scatter.zoom_to_table(), while_typing=False),
        Shortcut("V", "Table: only features in the visible region of the plot (on / off)", ("<Key-v>", "<Key-V>"),
                 lambda app: app.scatter.toggle_only_visible(), while_typing=False),
        Shortcut("M", "Table: only features marked in the plot (on / off)", ("<Key-m>", "<Key-M>"),
                 lambda app: app.scatter.toggle_only_marked(), while_typing=False),
    ]),
    ("Scatter plot, mouse over the plot", [Shortcut(keys, description) for _, keys, description, _ in PLOT_KEYS]),
    ("Annotation manager (Ctrl+M)", [
        Shortcut("Space / double-click", "Keep / discard the selected annotations"),
        Shortcut("Delete", "Discard the selected annotations"),
        Shortcut("Insert", "Keep the selected annotations"),
        Shortcut("Ctrl+A", "Select all"),
        Shortcut("Ctrl+C", "Copy the selected rows"),
        Shortcut("Ctrl+S", "Write the changes to the measured database"),
        Shortcut("Click a row", "Show its feature (table, overview, scatter plot)"),
    ]),
    ("Dialogs", [
        Shortcut("Enter", "Add (manual annotation, relation between features)"),
    ]),
]

# (where, mouse action, effect)
MOUSE_ACTIONS: list[tuple[str, str, str]] = [
    ("Scatter plot", "Click", "Select the feature: its best match in the table, otherwise the feature alone"),
    ("Scatter plot", "Shift+click", "Compare two features (the older one is replaced)"),
    ("Scatter plot", "Ctrl+click", "Add to / remove from a series (after releasing Ctrl: new series)"),
    ("Results table", "Click a row", "Show the match (the plot pans to the feature if needed)"),
    ("Results table", "Shift+click a row", "Add its feature to the comparison"),
    ("Results table", "Ctrl+click a row", "Add its feature to / remove it from a series (as in the scatter plot)"),
    ("Results table", "Click a header", "Sort: ascending → descending → default order"),
    ("Results table", "Ctrl+click a header", "Further sort key"),
    ("Results table", "Right click", "Copy the cell, the row or the table"),
    ("Annotation manager", "Click a header", "Sort, as in the results table (Ctrl+click: further sort key)"),
    ("Overview", "Click an MS2 peak", "Annotate the peak"),
    ("Everywhere", "Click a section header", "Collapse / expand the section"),
]

_TEXT_FIELDS = (tk.Entry, ttk.Entry, tk.Text, tk.Spinbox)


def bind_shortcuts(app: "AnnotationGUI"):
    """Bind the shortcuts to the window of the GUI."""
    window: tk.Misc = app.tk_root.winfo_toplevel()

    def run_shortcut(event, shortcut: Shortcut):
        if not app.tk_root.winfo_ismapped():  # e.g. another tab of the main window is shown
            return None
        typing: bool = isinstance(event.widget, _TEXT_FIELDS)
        if (not shortcut.while_typing and typing) or (shortcut.only_while_typing and not typing):
            return None
        shortcut.action(app)
        return "break"

    for _, shortcuts in SHORTCUTS:
        for shortcut in shortcuts:
            for sequence in shortcut.sequences:
                window.bind(sequence, lambda event, bound=shortcut: run_shortcut(event, bound), add="+")
    close_popups_with_escape(window)


def close_popups_with_escape(widget: tk.Misc):
    """
    Esc closes the pop-up window it is pressed in, as its close button would (e.g. the annotation manager still asks
    about unsaved changes). Applies to all windows of the application except the main window, bound once.
    """
    root: tk.Misc = widget._root()
    if getattr(root, "_escape_closes_popups", False):
        return
    root._escape_closes_popups = True

    def close(event):
        try:
            window = event.widget.winfo_toplevel()
        except (AttributeError, KeyError, tk.TclError):  # e.g. the window was closed by its own Esc binding already
            return None
        # the main window, and the list of an open combobox (its own toplevel, Esc closes only the list)
        if not isinstance(window, tk.Toplevel) or "popdown" in str(window):
            return None
        handler: str = window.protocol("WM_DELETE_WINDOW")
        if handler:
            window.tk.call(handler)
        else:
            window.destroy()
        return "break"

    root.bind_all("<Escape>", close, add="+")


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
    window.bind("<F1>", lambda event: window.destroy())  # Esc: see close_popups_with_escape

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
