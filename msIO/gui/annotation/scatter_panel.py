"""Scatter plot (m/z vs RT) of the measured features."""

import tkinter as tk
from tkinter import messagebox, ttk
from typing import TYPE_CHECKING

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from msIO.gui.annotation.common import (
    COMPARISON_COLORS, CollapsibleSection, find_feature_pairs, pack_canvas, parse_mass_differences, parse_masses
)

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI

# matplotlib key names of + and - (main keyboard and numpad) for zooming
ZOOM_IN_KEYS: set[str] = {"+", "plus", "add", "kp_add"}
ZOOM_OUT_KEYS: set[str] = {"-", "minus", "subtract", "kp_subtract"}
# colors of the lines connecting features with marked mass differences (one per entry)
MASS_DIFFERENCE_COLORS: list[str] = ["C3", "C4", "C5", "C6", "C8", "C9"]


class ScatterPanel(CollapsibleSection):
    """
    Scatter plot of all measured features with markers (annotated features, filtered matches, selected and compared
    features), highlighted features (MS2 fragments / neutral losses) and lines between features with marked mass
    differences. Click: select / show a feature, shift+click: compare (see AnnotationGUI.on_scatter_click),
    + / -: zoom.
    """

    ZOOM_FACTOR: float = 1.5

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, "Scatter plot m/z vs RT (click: select in table, shift-click two features: compare)",
                         on_toggle=lambda expanded: app.on_toggle_plot_section())
        self.app: "AnnotationGUI" = app
        self.tk_show_annotations: tk.BooleanVar = tk.BooleanVar(app.tk_root, value=True)
        self.tk_circle_matches: tk.BooleanVar = tk.BooleanVar(app.tk_root, value=True)
        # MS2 peaks below this intensity are ignored when searching fragments and losses
        self.tk_min_peak_intensity: tk.DoubleVar = tk.DoubleVar(app.tk_root, value=100)
        self.tk_rt_window: tk.DoubleVar = tk.DoubleVar(app.tk_root, value=0.2)

        self.canvas: FigureCanvasTkAgg = None
        self.toolbar: NavigationToolbar2Tk = None
        self.ax = None
        self.highlighted_f_ids: list[int] = None  # features with the requested fragments / losses
        # pairs of features with a marked mass difference: (entry, [(feature, partner), ...]) per entry
        self.mass_difference_pairs: list[tuple[str, list[tuple[int, int]]]] = []
        self.filtered_features_marker = None  # black circles around features with matches in the filtered table
        self.selected_feature_marker = None  # red circle around the selected feature
        self.compared_features_markers: list = []  # squares around the compared features
        self.mass_difference_lines: list = []  # lines and their end caps
        self.annotated_features_markers: list = []  # green / yellow points of annotated features
        self.focus_before_plot: tk.Misc = None  # gets the focus back when the mouse leaves the plot
        # annotated MS2 peaks of the measured database offered in the fragment list: {shown text: inserted value}
        self.fragment_choices: dict[str, str] = {}
        self.fragments_text_before_selection: str = ""

        tk_controls = ttk.Frame(self.body)
        tk_controls.pack(fill="x")
        ttk.Checkbutton(tk_controls, text="Show annotations", variable=self.tk_show_annotations,
                        command=self.redraw).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(tk_controls, text="Circle matches", variable=self.tk_circle_matches,
                        command=self.mark_filtered_features).grid(row=0, column=1, sticky="w", padx=(4, 10))
        ttk.Label(tk_controls, text="Min. peak intensity").grid(row=0, column=2, sticky="w")
        ttk.Entry(tk_controls, textvariable=self.tk_min_peak_intensity, width=10).grid(row=0, column=3, sticky="w",
                                                                                        padx=(4, 10))
        ttk.Label(tk_controls, foreground="gray", text="(fragments and losses; + / - over the plot: zoom)").grid(
            row=0, column=4, columnspan=4, sticky="w")

        # highlight features with MS2 fragments / neutral losses, entries are m/z values or formulas
        ttk.Label(tk_controls, text="Fragments (a; b)").grid(row=1, column=0, sticky="w", pady=2)
        # can be typed in, the list holds the annotated peaks of the measured database (selecting one appends it)
        self.tk_fragments_entry = ttk.Combobox(tk_controls, width=22, postcommand=self._remember_fragments_text)
        self.tk_fragments_entry.grid(row=1, column=1, sticky="w", padx=(4, 10))
        self.tk_fragments_entry.bind("<<ComboboxSelected>>", self.on_select_fragment)
        ttk.Label(tk_controls, text="Losses (a; b)").grid(row=1, column=2, sticky="w")
        self.tk_losses_entry = ttk.Entry(tk_controls, width=22)
        self.tk_losses_entry.grid(row=1, column=3, sticky="w", padx=(4, 10))
        ttk.Label(tk_controls, text="Match").grid(row=1, column=4, sticky="w")
        self.tk_highlight_mode: ttk.Combobox = ttk.Combobox(tk_controls, state="readonly", values=["all", "any"],
                                                            width=4)
        self.tk_highlight_mode.set("all")
        self.tk_highlight_mode.grid(row=1, column=5, sticky="w", padx=(4, 10))
        ttk.Button(tk_controls, text="Highlight", command=self.on_click_highlight).grid(row=1, column=6, padx=3)
        ttk.Button(tk_controls, text="Remove highlight", command=self.remove_highlight).grid(row=1, column=7)
        for tk_entry in (self.tk_fragments_entry, self.tk_losses_entry):
            tk_entry.bind("<Return>", lambda event: self.on_click_highlight())

        # connect pairs of features with mass differences or adduct exchanges, matched with the precursor tolerance
        ttk.Label(tk_controls, text="Mass differences").grid(row=2, column=0, sticky="w", pady=2)
        self.tk_mass_differences_entry = ttk.Entry(tk_controls)
        self.tk_mass_differences_entry.grid(row=2, column=1, columnspan=3, sticky="ew", padx=(4, 10))
        self.tk_mass_differences_entry.bind("<Return>", lambda event: self.on_click_mark_mass_differences())
        ttk.Label(tk_controls, text="RT window (min)").grid(row=2, column=4, sticky="w")
        ttk.Entry(tk_controls, textvariable=self.tk_rt_window, width=5).grid(row=2, column=5, sticky="w", padx=(4, 10))
        ttk.Button(tk_controls, text="Mark", command=self.on_click_mark_mass_differences).grid(row=2, column=6, padx=3)
        ttk.Button(tk_controls, text="Remove marks", command=self.remove_mass_differences).grid(row=2, column=7)
        ttk.Label(tk_controls, foreground="gray",
                  text="e.g. NH3; C[13] - C[12]; [M+H]+>[M+NH4]+ (adduct exchange); fragments are ions, "
                       "e.g. C5H15NO4P+").grid(row=3, column=0, columnspan=8, sticky="w")

        self.tk_plot_frame = ttk.Frame(self.body)
        self.tk_plot_frame.pack(fill="both", expand=True)

    # ---- drawing ----------------------------------------------------------------------------------------------------
    def reset(self):
        """New measured database: no highlights and marks, full view."""
        self.highlighted_f_ids = None
        self.mass_difference_pairs = []
        self.draw(keep_view=False)
        self.update_fragment_choices()

    def _destroy_plot(self):
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        if self.toolbar is not None:
            self.toolbar.destroy()
            self.toolbar = None
        self.ax = None
        self.filtered_features_marker = None
        self.selected_feature_marker = None
        self.compared_features_markers = []
        self.mass_difference_lines = []
        self.annotated_features_markers = []

    def draw(self, keep_view: bool = True):
        """
        (Re)draw the scatter plot of the measurement including annotations, highlights and markers.
        With keep_view, the zoomed / panned region of the previous plot is kept (not for a new database).
        """
        view = None
        if keep_view and self.ax is not None:
            view = self.ax.get_xlim(), self.ax.get_ylim()
        self._destroy_plot()

        fig = Figure(figsize=(7, 6))
        self.ax = fig.add_subplot(111)
        try:
            self.app.meas_db.plot_scatter(add_annotations=self.tk_show_annotations.get(),
                                          highlight_features=self.highlighted_f_ids, ax=self.ax)
        except Exception as error:
            self.ax = None
            messagebox.showerror("Error", f"Could not render the scatter plot:\n{error}")
            return
        self.ax.set_xlabel("Retention time (min)")
        self.ax.set_ylabel("m/z")
        fig.tight_layout()

        self.canvas = FigureCanvasTkAgg(fig, master=self.tk_plot_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.tk_plot_frame)
        self.toolbar.update()
        self.toolbar.pack()
        pack_canvas(self.canvas)
        self.canvas.mpl_connect('button_press_event', self.on_click_plot)
        self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.mpl_connect('figure_enter_event', self.on_enter)
        self.canvas.mpl_connect('figure_leave_event', self.on_leave)

        self.mark_annotated_features()
        self.mark_filtered_features()
        if (match := self.app.results.selected_match()) is not None:
            self.mark_selected_feature(int(match.feature_id_meas))
        self.mark_compared_features()
        self.mark_mass_differences()
        if view is not None:
            # the full view first, so that the home button of the toolbar still shows all features
            self.toolbar.push_current()
            self.ax.set_xlim(*view[0])
            self.ax.set_ylim(*view[1])
            self.toolbar.push_current()
        self.canvas.draw()

    def redraw(self):
        """Redraw the scatter plot, if there is one (i.e. a measured database is loaded)."""
        if self.ax is not None:
            self.draw()

    def _update_legend(self):
        # only the markers have labels, so the legend lists exactly the markers that are shown
        if self.ax.get_legend_handles_labels()[0]:
            self.ax.legend(loc='upper right')
        elif (legend := self.ax.get_legend()) is not None:
            legend.remove()
        self.canvas.draw_idle()

    def _points(self, f_ids) -> list[tuple[float, float]]:
        return [p for f_id in f_ids if (p := self.app.get_rt_mz(f_id)) is not None]

    # ---- markers ----------------------------------------------------------------------------------------------------
    def mark_annotated_features(self):
        """
        Color features with manual annotations (in the database or to add) green, highlighted ones in yellow
        (red + green), independent of whether the annotation labels are shown.
        """
        if self.ax is None:
            return
        for artist in self.annotated_features_markers:
            artist.remove()
        self.annotated_features_markers = []
        annotated: set[int] = self.app.annotations.annotated_feature_ids()
        highlighted: set[int] = set(self.highlighted_f_ids or [])
        for f_ids, color, label in ((annotated - highlighted, "green", "Annotated"),
                                    (annotated & highlighted, "gold", "Annotated + highlighted")):
            if points := self._points(f_ids):
                rts, mzs = zip(*points)
                # same size as the highlighted features of FeatureManagerDB.plot_scatter
                self.annotated_features_markers.append(self.ax.scatter(rts, mzs, c=color, s=20, zorder=3, label=label))
        self._update_legend()

    def mark_filtered_features(self):
        """Circle features of the filtered table in black (only if a filter removes rows and circles are enabled)."""
        if self.ax is None:
            return
        if self.filtered_features_marker is not None:
            self.filtered_features_marker.remove()
        self.filtered_features_marker = None
        if self.tk_circle_matches.get() and (points := self._points(self.app.get_filtered_feature_ids())):
            rts, mzs = zip(*points)
            self.filtered_features_marker = self.ax.scatter(rts, mzs, s=120, facecolors='none', edgecolors='black',
                                                            linewidths=2, zorder=10, label='Library match')
        self._update_legend()

    def mark_selected_feature(self, f_id_meas: int):
        """Circle the given feature in red."""
        if self.ax is None or (point := self.app.get_rt_mz(f_id_meas)) is None:
            return
        if self.selected_feature_marker is not None:
            self.selected_feature_marker.remove()
        self.selected_feature_marker = self.ax.scatter([point[0]], [point[1]], s=300, facecolors="none",
                                                       edgecolors="red", linewidths=2, zorder=10,
                                                       label='Selected feature')
        self._update_legend()

    def mark_compared_features(self):
        """Put squares around the compared features (AnnotationGUI.compared_f_ids)."""
        if self.ax is None:
            return
        for artist in self.compared_features_markers:
            artist.remove()
        self.compared_features_markers = []
        for f_id, color in zip(self.app.compared_f_ids, COMPARISON_COLORS):
            if (point := self.app.get_rt_mz(f_id)) is None:
                continue
            self.compared_features_markers.append(self.ax.scatter(
                [point[0]], [point[1]], s=200, marker="s", facecolors="none", edgecolors=color, linewidths=2,
                zorder=11, label=f"Compared feature {f_id}"
            ))
        self._update_legend()

    def mark_mass_differences(self):
        """Draw a line between the features of each pair with a marked mass difference."""
        if self.ax is None:
            return
        for artist in self.mass_difference_lines:
            artist.remove()
        self.mass_difference_lines = []
        for i, (entry, pairs) in enumerate(self.mass_difference_pairs):
            color: str = MASS_DIFFERENCE_COLORS[i % len(MASS_DIFFERENCE_COLORS)]
            segments = [(self.app.get_rt_mz(f_id_a), self.app.get_rt_mz(f_id_b)) for f_id_a, f_id_b in pairs]
            lines = LineCollection(segments, colors=color, linewidths=2.5, zorder=12, label=f"{entry} ({len(pairs)})")
            self.ax.add_collection(lines)
            self.mass_difference_lines.append(lines)
            # mass differences are small compared to the m/z range, the end caps keep short lines visible
            if segments:
                ends = np.array([point for segment in segments for point in segment])
                self.mass_difference_lines.append(self.ax.scatter(ends[:, 0], ends[:, 1], marker="_", s=120,
                                                                  linewidths=2.5, color=color, zorder=12))
        self._update_legend()

    # ---- highlighted features ---------------------------------------------------------------------------------------
    def on_click_highlight(self):
        """
        Highlight features whose MS2 spectrum contains the entered fragments and neutral losses.

        Peaks are matched with the MS2 fragment tolerance, peaks below the minimum peak intensity are ignored.
        With "all", a feature needs every entered fragment and loss, with "any" one of them is enough.
        Empty fields remove the highlight.
        """
        if self.app.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        try:
            fragments = parse_masses(self.tk_fragments_entry.get())
            losses = parse_masses(self.tk_losses_entry.get())
            mz_tol_da: float = self.app.matching.ms2_tolerance_da()
            min_intensity: float = self.tk_min_peak_intensity.get()
        except ValueError as error:
            messagebox.showerror("Invalid input", str(error))
            return
        except tk.TclError:
            messagebox.showerror("Invalid input", "Check the MS2 fragment tolerance and the minimum peak intensity.")
            return
        if not (fragments or losses):
            self.remove_highlight()
            return

        meas_db = self.app.meas_db
        try:
            f_id_sets: list[set[int]] = (
                    [set(meas_db.filter_by_fragment_mz(mass, mz_tol_da, min_intensity)) for _, mass in fragments]
                    + [set(meas_db.filter_by_neutral_loss(mass, mz_tol_da, min_intensity)) for _, mass in losses]
            )
        except Exception as error:
            messagebox.showerror("Error", f"Could not search for fragments and losses:\n{error}")
            return
        mode: str = self.tk_highlight_mode.get()
        combined = set.intersection(*f_id_sets) if mode == "all" else set.union(*f_id_sets)
        self.highlighted_f_ids = sorted(combined)
        self.draw()

        def describe(entries: list[tuple[str, float]]) -> str:
            return ", ".join(entry if entry == f"{mass:g}" else f"{entry} ({mass:.4f})" for entry, mass in entries)

        searched = "; ".join(part for part in (
            f"fragments {describe(fragments)}" if fragments else "",
            f"losses {describe(losses)}" if losses else "",
        ) if part)
        self.app.set_status(f"Highlighted {len(self.highlighted_f_ids)} features with {mode} of {searched} "
                            f"(peaks above {min_intensity:g}).")

    def remove_highlight(self):
        self.highlighted_f_ids = None
        self.redraw()

    def update_fragment_choices(self):
        """Fill the fragment list with the annotated MS2 peaks of the measured database, most frequent first."""
        self.fragment_choices = {}
        try:
            annotations = self.app.meas_db.get_peak_annotations(level=2) if self.app.meas_db is not None else []
        except Exception:  # e.g. not an msIO database
            annotations = []
        for annotation, mz, n_features in annotations:
            try:  # formulas (e.g. C5H15NO4P+) are inserted as they are, other names as their m/z
                parse_masses(annotation)
                value = annotation
            except ValueError:
                value = f"{mz:.4f}"
            label = f"{annotation}  (m/z {mz:.4f}, {n_features} feature{'s' if n_features > 1 else ''})"
            self.fragment_choices[label] = value
        self.tk_fragments_entry.config(values=list(self.fragment_choices))

    def _remember_fragments_text(self):
        # the combobox replaces the text with the selected entry, keep what was typed before
        self.fragments_text_before_selection = self.tk_fragments_entry.get()

    def on_select_fragment(self, event=None):
        """Append the selected annotated fragment to the fragments typed before."""
        value: str = self.fragment_choices.get(self.tk_fragments_entry.get(), self.tk_fragments_entry.get())
        before: str = self.fragments_text_before_selection.strip().rstrip(";").strip()
        self.tk_fragments_entry.set(f"{before}; {value}" if before else value)
        self.tk_fragments_entry.icursor("end")

    # ---- mass differences -------------------------------------------------------------------------------------------
    def on_click_mark_mass_differences(self):
        """
        Connect pairs of features whose m/z differ by the entered mass differences or adduct exchanges.

        m/z values are matched with the precursor tolerance (stricter of ppm and mDa), the retention times of
        both features have to be within the RT window. Empty field removes the marks.
        """
        if self.app.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        try:
            differences = parse_mass_differences(self.tk_mass_differences_entry.get())
            max_dmz_ppm, max_dmz_da = self.app.matching.precursor_tolerances()
            rt_window: float = self.tk_rt_window.get()
        except ValueError as error:
            messagebox.showerror("Invalid input", str(error))
            return
        except tk.TclError:
            messagebox.showerror("Invalid input", "Check the precursor tolerances and the RT window.")
            return
        if not differences:
            self.remove_mass_differences()
            return

        points = {f_id: p for f_id in self.app.meas_db.feature_ids if (p := self.app.get_rt_mz(f_id)) is not None}
        f_ids = np.array(list(points.keys()))
        rts, mzs = np.array(list(points.values())).T
        self.mass_difference_pairs = [
            (entry, find_feature_pairs(f_ids, rts, mzs, expected_partner_mz, max_dmz_ppm, max_dmz_da, rt_window))
            for entry, expected_partner_mz in differences
        ]
        self.mark_mass_differences()
        counts = ", ".join(f"{entry}: {len(pairs)}" for entry, pairs in self.mass_difference_pairs)
        self.app.set_status(f"Pairs of features within {rt_window:g} min: {counts}.")

    def remove_mass_differences(self):
        self.mass_difference_pairs = []
        self.mark_mass_differences()

    def find_marked_mass_differences(self, f_id_a: int, f_id_b: int) -> list[str]:
        """Marked mass differences that connect two features (in either direction)."""
        return [entry for entry, pairs in self.mass_difference_pairs
                if (f_id_a, f_id_b) in pairs or (f_id_b, f_id_a) in pairs]

    # ---- mouse and keys ---------------------------------------------------------------------------------------------
    def find_closest_feature(self, x: float, y: float) -> int | None:
        """
        Feature closest to a point of the plot (highlighted or not).
        RT and m/z have different scales, so the distance is measured relative to the visible axis range.
        """
        points = {f_id: p for f_id in self.app.meas_db.feature_ids if (p := self.app.get_rt_mz(f_id)) is not None}
        if not points:
            return None
        f_ids = np.array(list(points.keys()))
        rts, mzs = np.array(list(points.values())).T
        x_min, x_max = self.ax.get_xlim()
        y_min, y_max = self.ax.get_ylim()
        distances = ((rts - x) / (x_max - x_min)) ** 2 + ((mzs - y) / (y_max - y_min)) ** 2
        return int(f_ids[np.argmin(distances)])

    def on_click_plot(self, event):
        # ignore clicks outside the axes or while zooming / panning with the toolbar
        if event.inaxes is not self.ax or self.toolbar.mode:
            return
        if (closest_f_id := self.find_closest_feature(event.xdata, event.ydata)) is None:
            return
        # modifiers come from the state of this click; event.key is the last key seen by the canvas and can get
        # stuck at "shift" when the key is released while another widget has the focus
        self.app.on_scatter_click(closest_f_id, shift="shift" in event.modifiers)

    def on_enter(self, event):
        """
        Keys only reach the canvas if it has the focus: take it while the mouse is over the plot and give it back to
        the previous widget (e.g. a filter field) when the mouse leaves.
        """
        try:
            focused = self.app.tk_root.focus_get()
        except KeyError:  # tkinter can't resolve some internal widgets, e.g. the open list of a combobox
            focused = None
        canvas_widget = self.canvas.get_tk_widget()
        if focused is not canvas_widget:
            self.focus_before_plot = focused
            canvas_widget.focus_set()

    def on_leave(self, event):
        canvas_widget = self.canvas.get_tk_widget()
        if self.focus_before_plot is not None and self.app.tk_root.focus_get() is canvas_widget:
            try:
                self.focus_before_plot.focus_set()
            except tk.TclError:  # widget was destroyed in the meantime
                pass
        self.focus_before_plot = None

    def on_key(self, event):
        """+ / - zoom in / out by ZOOM_FACTOR, centered on the mouse position."""
        # key names depend on the keyboard and on whether the key produces a character (e.g. numpad)
        if event.inaxes is not self.ax or event.key not in ZOOM_IN_KEYS | ZOOM_OUT_KEYS:
            return
        scale: float = 1 / self.ZOOM_FACTOR if event.key in ZOOM_IN_KEYS else self.ZOOM_FACTOR
        # the point under the mouse stays where it is
        x_min, x_max = self.ax.get_xlim()
        y_min, y_max = self.ax.get_ylim()
        self.ax.set_xlim(event.xdata - (event.xdata - x_min) * scale, event.xdata + (x_max - event.xdata) * scale)
        self.ax.set_ylim(event.ydata - (event.ydata - y_min) * scale, event.ydata + (y_max - event.ydata) * scale)
        # makes the zoom steps available to the back / forward buttons of the toolbar
        self.toolbar.push_current()
        self.canvas.draw_idle()
