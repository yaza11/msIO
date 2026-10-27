"""Overview of a library match, a single feature or two compared features (figure: msIO.gui.feature_overview)."""

import re
import threading
import tkinter as tk
from tkinter import messagebox, simpledialog, ttk
from typing import TYPE_CHECKING

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from msIO.feature_managers.db import FeatureManagerDB, Library
from msIO.gui.annotation.common import CollapsibleSection, pack_canvas
from msIO.gui.feature_overview import (
    OVERVIEW_PANELS, OverviewData, OverviewSettings, PlottedSpectrum, SpectrumKey, build_overview_figure,
    gather_feature, gather_library, peak_key
)

if TYPE_CHECKING:
    from msIO.gui.annotation.app import AnnotationGUI

# labels of the checkboxes that show / hide the plots of the overview
PANEL_LABELS: dict[str, str] = {"ms2": "MS2", "ms1": "MS1", "intensities": "Intensities", "structure": "Structure"}


class OverviewPanel(CollapsibleSection):
    """
    Shows an overview: "match" (measured feature and library match), "feature" (one measured feature) or
    "comparison" (two measured features). Gathering the data and building the figure run in a worker thread.

    The controls set the shown feature(s) by ID, hide the isotope peaks of the precursor, set the cutoff of the peak
    annotations and select the plots. Clicking an MS2 peak edits its annotation.
    """

    DEFAULT_TITLE: str = "Overview (select a match in the table or features in the scatter plot)"

    def __init__(self, parent: tk.Misc, app: "AnnotationGUI"):
        super().__init__(parent, self.DEFAULT_TITLE, on_toggle=lambda expanded: app.on_toggle_plot_section())
        self.app: "AnnotationGUI" = app
        self.tk_hide_isotope_peaks: tk.BooleanVar = tk.BooleanVar(app.tk_root, value=False)
        self.tk_annotation_cutoff: tk.DoubleVar = tk.DoubleVar(app.tk_root, value=0.3)
        self.tk_show_panels: dict[str, tk.BooleanVar] = {
            panel: tk.BooleanVar(app.tk_root, value=True) for panel in OVERVIEW_PANELS
        }
        self.annotation_cutoff_job: str = None

        self.canvas: FigureCanvasTkAgg = None
        self.toolbar: NavigationToolbar2Tk = None
        self.running: bool = False
        # latest request: (request number, key, match result, feature annotations, data to reuse, settings)
        self.pending: tuple = None
        # incremented for every new content of the panel, so a slow overview can't replace newer content
        self.request_number: int = 0
        # (mode, feature ids, library feature id) of the shown overview, requesting it again doesn't replot it
        self.shown_key: tuple = None
        self.data: OverviewData = None  # data of the shown overview, reused when only settings change
        self.spectra: list[PlottedSpectrum] = None  # plotted MS2 spectra, for clicks on peaks
        # annotations of library MS2 peaks set in this session (the library is not changed):
        # {("lib", feature id): {m/z (4 decimals): annotation}}
        self.peak_annotations: dict[SpectrumKey, dict[float, str]] = {}

        tk_controls = ttk.Frame(self.body)
        tk_controls.pack(fill="x")
        ttk.Label(tk_controls, text="Feature ID(s)").pack(side="left")
        self.tk_feature_ids_entry = ttk.Entry(tk_controls, width=14)
        self.tk_feature_ids_entry.pack(side="left", padx=(4, 0))
        self.tk_feature_ids_entry.bind("<Return>", lambda event: self.on_submit_feature_ids())
        ttk.Button(tk_controls, text="Show", command=self.on_submit_feature_ids).pack(side="left", padx=(3, 15))
        ttk.Checkbutton(tk_controls, text="Hide precursor isotope peaks (±4 Da)", variable=self.tk_hide_isotope_peaks,
                        command=self.redraw).pack(side="left")
        ttk.Label(tk_controls, text="Annotation cutoff").pack(side="left", padx=(15, 3))
        ttk.Scale(tk_controls, from_=0, to=1, variable=self.tk_annotation_cutoff, length=120,
                  command=self.on_change_annotation_cutoff).pack(side="left")
        self.tk_annotation_cutoff_label = ttk.Label(tk_controls, text=f"{self.tk_annotation_cutoff.get():.2f}", width=4)
        self.tk_annotation_cutoff_label.pack(side="left", padx=(3, 0))

        tk_panel_controls = ttk.Frame(self.body)
        tk_panel_controls.pack(fill="x", pady=(2, 0))
        ttk.Label(tk_panel_controls, text="Plots:").pack(side="left")
        for panel, variable in self.tk_show_panels.items():
            ttk.Checkbutton(tk_panel_controls, text=PANEL_LABELS[panel], variable=variable,
                            command=self.redraw).pack(side="left", padx=(6, 0))
        ttk.Label(tk_panel_controls, foreground="gray", text="(click an MS2 peak to annotate it)").pack(
            side="left", padx=15)

        self.tk_plot_frame = ttk.Frame(self.body)
        self.tk_plot_frame.pack(fill="both", expand=True)

    # ---- requests ---------------------------------------------------------------------------------------------------
    def _get_settings(self) -> OverviewSettings | None:
        try:
            ms2_tolerance_da: float = self.app.matching.ms2_tolerance_da()
            ms1_tolerance_da: float = self.app.matching.precursor_tolerances()[1]
            annotation_cutoff: float = self.tk_annotation_cutoff.get()
        except tk.TclError:
            messagebox.showerror("Invalid Tolerance", "Check the m/z tolerances.")
            return None
        return OverviewSettings(
            hide_isotope_peaks=self.tk_hide_isotope_peaks.get(),
            annotation_cutoff=annotation_cutoff,
            ms2_tolerance_da=ms2_tolerance_da,
            ms1_tolerance_da=ms1_tolerance_da,
            # copy, the worker reads it while the user may edit annotations
            peak_annotations={key: dict(annotations) for key, annotations in self.peak_annotations.items()},
            panels=tuple(panel for panel, variable in self.tk_show_panels.items() if variable.get()),
        )

    def request(self, mode: str, f_ids: list[int], f_id_lib: int = None, reuse_data: bool = False):
        """
        Show an overview: "match" (feature f_ids[0] and library feature f_id_lib), "feature" (f_ids[0]) or
        "comparison" (both f_ids).

        Only one worker runs at a time; if another overview is requested meanwhile, only the latest one is shown
        afterwards. With reuse_data, only the settings changed (cutoff, isotope filter, peak annotations, plots) and
        the data of the shown overview is used again.
        """
        if (settings := self._get_settings()) is None:
            return
        key: tuple = (mode, tuple(f_ids), f_id_lib)
        data: OverviewData | None = self.data if reuse_data and key == self.shown_key else None
        match_result = None
        if mode == "match":
            match_result = next(m for m in self.app.matches[f_ids[0]] if m['feature_id'] == f_id_lib)
        annotations = {f_id: self.app.annotations.feature_annotation(f_id) for f_id in f_ids}

        self.request_number += 1
        self.shown_key = key
        self.pending = (self.request_number, key, match_result, annotations, data, settings)
        if not reuse_data:
            self.app.set_status(f"Plotting {mode} of feature(s) {', '.join(map(str, f_ids))} ...")
        if not self.running:
            self._start_worker()

    def redraw(self):
        """Redraw the shown overview with the current settings."""
        if self.shown_key is not None and self.data is not None:
            self.request(*self.shown_key, reuse_data=True)

    def clear(self):
        """Remove the overview (e.g. new database), an overview that is still being built is dropped."""
        self.request_number += 1
        self.pending = None
        self.shown_key = None
        self.data = self.spectra = None
        self._destroy_plot()
        self.set_title(self.DEFAULT_TITLE)

    def _start_worker(self):
        request_number, key, match_result, annotations, data, settings = self.pending
        self.pending = None
        self.running = True
        threading.Thread(
            target=self._worker,
            args=(self.app.meas_db, self.app.lib_db, self.app.ms2_spectra, request_number, key, match_result,
                  annotations, data, settings),
            daemon=True
        ).start()

    def _worker(self, meas_db: FeatureManagerDB, lib_db: Library, ms2_spectra: dict | None, request_number: int,
                key: tuple, match_result: dict | None, annotations: dict, data: OverviewData | None,
                settings: OverviewSettings):
        """Worker thread: gather the data (unless given) and build the figure (not attached to tk yet)."""
        fig, spectra, error_text = None, None, None
        try:
            if data is None:
                mode, f_ids, f_id_lib = key
                data = OverviewData(
                    mode=mode,
                    features=[gather_feature(meas_db, f_id, ms2_spectra, annotations[f_id]) for f_id in f_ids],
                    library=gather_library(lib_db, f_id_lib) if mode == "match" else None,
                    match_result=match_result,
                )
            fig = Figure(figsize=(12, 10))
            spectra = build_overview_figure(fig, data, settings)
        except Exception as error:
            error_text = f"{type(error).__name__}: {error}"
        self.app.call_in_main_thread(
            lambda: self._finish(request_number, key, data, fig, spectra, settings, error_text))

    def _finish(self, request_number: int, key: tuple, data: OverviewData | None, fig: Figure | None,
                spectra: list[PlottedSpectrum] | None, settings: OverviewSettings, error_text: str | None):
        """Main thread: show the figure built by the worker, unless newer content was requested in the meantime."""
        self.running = False
        if self.pending is not None:
            self._start_worker()
            return
        if request_number != self.request_number:
            return

        if error_text is not None:
            self.shown_key = None  # requesting it again retries
            self.app.set_status("Could not plot the overview.")
            messagebox.showerror("Plot Error", f"Could not plot the overview:\n{error_text}")
            return

        self.data, self.spectra = data, spectra
        mode, f_ids, f_id_lib = key
        title = {
            "match": f"Match: feature {f_ids[0]} -> library feature {f_id_lib}",
            "feature": f"Feature {f_ids[0]}",
            "comparison": f"Comparison: feature {f_ids[0]} vs {f_ids[-1]}",
        }[mode]
        self._show(fig, title)
        status: str = f"Plotted {title[0].lower()}{title[1:]}."
        if settings.hide_isotope_peaks and "ms2" in settings.panels:
            status += f" Hid {self._count_hidden_isotope_peaks(data, spectra)} MS2 peak(s) within ±4 Da of the precursor."
        self.app.set_status(status)

    @staticmethod
    def _count_hidden_isotope_peaks(data: OverviewData, spectra: list[PlottedSpectrum]) -> int:
        """Number of MS2 peaks removed by the isotope filter (all shown spectra)."""
        n_peaks: dict[SpectrumKey, int] = {("meas", f.f_id): len(f.ms2.mzs) for f in data.features if f.ms2 is not None}
        if data.library is not None and data.library.ms2 is not None:
            n_peaks[("lib", data.library.f_id)] = len(data.library.ms2.mzs)
        return sum(n_peaks[spectrum.key] - len(spectrum.mzs) for spectrum in spectra if spectrum.key in n_peaks)

    # ---- canvas -----------------------------------------------------------------------------------------------------
    def _show(self, fig: Figure, title: str):
        """Replace the content of the panel."""
        self._destroy_plot()
        self.set_title(title)
        self.canvas = FigureCanvasTkAgg(fig, master=self.tk_plot_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.tk_plot_frame)
        self.toolbar.update()
        self.toolbar.pack()
        pack_canvas(self.canvas)
        self.canvas.mpl_connect('button_press_event', self.on_click_plot)
        self.canvas.draw()

    def _destroy_plot(self):
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        if self.toolbar is not None:
            self.toolbar.destroy()
            self.toolbar = None

    # ---- controls ---------------------------------------------------------------------------------------------------
    def on_change_annotation_cutoff(self, value: str):
        """Slider moved: update the label right away, redraw once the slider rests."""
        self.tk_annotation_cutoff_label.config(text=f"{float(value):.2f}")
        if self.annotation_cutoff_job is not None:
            self.app.tk_root.after_cancel(self.annotation_cutoff_job)
        self.annotation_cutoff_job = self.app.tk_root.after(300, self._apply_annotation_cutoff)

    def _apply_annotation_cutoff(self):
        self.annotation_cutoff_job = None
        self.redraw()

    def on_submit_feature_ids(self):
        """Feature ID field: one ID shows that feature, two IDs (separated by , ; or space) compare them."""
        if self.app.meas_db is None:
            messagebox.showerror("Error", "Please load a measured database first.")
            return
        try:
            f_ids = [int(part) for part in re.split(r"[,;\s]+", self.tk_feature_ids_entry.get().strip()) if part]
        except ValueError:
            messagebox.showerror("Invalid input", "Feature IDs must be whole numbers.")
            return
        if not 1 <= len(f_ids) <= 2 or len(set(f_ids)) != len(f_ids):
            messagebox.showerror("Invalid input", "Enter one feature ID or two different ones.")
            return
        if missing := [f_id for f_id in f_ids if f_id not in self.app.meas_db.mzs]:
            messagebox.showerror("Invalid input", f"There is no measured feature {', '.join(map(str, missing))}.")
            return
        self.app.show_feature_ids(f_ids)

    # ---- peak annotations -------------------------------------------------------------------------------------------
    def on_click_plot(self, event):
        """Click on an MS2 peak: edit its annotation (empty annotation shows the m/z again)."""
        if self.toolbar.mode or event.button != 1 or event.inaxes is None or self.spectra is None:
            return
        # upper or lower (mirror) spectrum, depending on where the click is
        candidates = [spectrum for spectrum in self.spectra
                      if spectrum.ax is event.inaxes and (spectrum.heights[0] > 0) == (event.ydata >= 0)]
        if not candidates:
            return
        spectrum = candidates[0]
        # closest peak in screen coordinates that reaches up to the click (small neighbors of a large peak are often
        # only a pixel away), at most a few pixels away
        x_pixels = event.inaxes.transData.transform(np.column_stack([spectrum.mzs, np.zeros(len(spectrum.mzs))]))[:, 0]
        distances = np.abs(x_pixels - event.x)
        distances[np.abs(spectrum.heights) + 0.1 < abs(event.ydata)] = np.inf
        idx: int = int(np.argmin(distances))
        if distances[idx] > 6:
            return

        mz: float = float(spectrum.mzs[idx])
        source, f_id = spectrum.key
        if source == "meas":
            # measured spectra: the annotation is stored in the measured database (and in the loaded spectra)
            peaks = self._get_loaded_peaks(f_id, mz)
            current: str = next((peak.annotation for peak in peaks if peak.annotation), "")
            where = "stored in the measured database"
        else:
            # library spectra are not changed, the annotation is only kept in this session
            current = self.peak_annotations.get(spectrum.key, {}).get(peak_key(mz), "")
            where = "only for this session"
        text = simpledialog.askstring(
            "Peak annotation", f"Annotation of the peak at m/z {mz:.4f} "
                               f"({'library feature' if source == 'lib' else 'feature'} {f_id}, {where}),\n"
                               f"empty to show the m/z again:",
            initialvalue=current, parent=self.app.tk_root
        )
        if text is None:  # cancelled
            return
        annotation: str | None = text.strip() or None
        if source == "meas":
            try:
                self.app.meas_db.set_peak_annotation(f_id, 2, mz, annotation)
            except Exception as error:
                messagebox.showerror("Peak annotation", f"Could not store the annotation:\n{error}")
                return
            for peak in peaks:
                peak.annotation = annotation
            self.app.scatter.update_fragment_choices()
        elif annotation:
            self.peak_annotations.setdefault(spectrum.key, {})[peak_key(mz)] = annotation
        else:
            self.peak_annotations.get(spectrum.key, {}).pop(peak_key(mz), None)
        self.redraw()

    def _get_loaded_peaks(self, f_id: int, mz: float) -> list:
        """Peaks at mz of the MS2 spectra of a measured feature that are loaded (overview and cache of all spectra)."""
        spectra = [feature.ms2 for feature in (self.data.features if self.data else []) if feature.f_id == f_id]
        if self.app.ms2_spectra is not None:
            spectra.append(self.app.ms2_spectra.get(f_id))
        peaks: dict[int, object] = {}  # by id, the overview usually uses the spectrum of the cache
        for spectrum in spectra:
            for peak in (spectrum.peaks if spectrum is not None else []):
                if peak_key(peak.mz) == peak_key(mz):
                    peaks[id(peak)] = peak
        return list(peaks.values())
