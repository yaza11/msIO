import functools
import sqlite3
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import messagebox, filedialog, ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from msIO.feature_managers.db import Library, FeatureManagerDB


def add_annotation_to_measured_db(f_id_meas: int, f_id_lib: int, results: pd.DataFrame, lib: Library, meas: FeatureManagerDB) -> None:
    """Add a matched library entry to the measured database."""
    ms2_score: float = results.loc[
        (results.feature_id_meas == f_id_meas) & (results.feature_id_lib == f_id_lib), 'ms2_score'
    ].squeeze()
    meas.add_annotations_from_library(lib, f_id_meas, f_id_lib, ms2_score)



class AnnotationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Annotation GUI")
        self.root.geometry("1500x750")
        self.meas_path = tk.StringVar()
        self.lib_path = tk.StringVar()
        self.f_id = tk.IntVar(value=1)
        self.max_dmz_ppm = tk.StringVar(value="10.0")
        self.max_dmz_mda = tk.StringVar(value="")
        self.max_dmz_ms2_mda = tk.StringVar(value="10")
        self.require_ms2 = tk.BooleanVar(value=False)

        self.meas = None
        self.lib = None
        self.current_matches = None
        self.canvas = None
        self.toolbar = None

        self.create_gui()

    def create_gui(self):

        #Seelecting database
        database_frame = ttk.LabelFrame(self.root, text="Database", padding=10)
        database_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(database_frame, text="Measured Database").grid(row=0, column=0, sticky="w")
        ttk.Entry(database_frame, textvariable=self.meas_path, width=55).grid(row=0, column=1, padx=5)
        ttk.Button(database_frame, text="Browse", command=self.load_meas_db).grid(row=0, column=2)

        ttk.Label(database_frame, text="Library Database").grid(row=1, column=0, sticky="w")
        ttk.Entry(database_frame, textvariable=self.lib_path, width=55).grid(row=1, column=1, padx=5)
        ttk.Button(database_frame, text="Browse", command=self.load_lib_db).grid(row=1, column=2)

        #Matching Parmeters
        #Parameter Frame
        parameter_frame = ttk.LabelFrame(self.root, text="Matching Parameters", padding=10)
        parameter_frame.pack(fill="x", padx=10, pady=5)

        # ttk.Label(parameter_frame, text="Feature ID (f_id):").grid(row=0, column=0, sticky="w")
        # ttk.Entry(parameter_frame, textvariable=self.f_id, width=10).grid(row=0, column=1, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="m/z tol. in ppm:").grid(row=0, column=0, sticky="w", padx=(10,0))
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_ppm, width=10).grid(row=0, column=1, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="m/z tol. in mDa:").grid(row=0, column=2, sticky="w", padx=(10,0))
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_mda, width=10).grid(row=0, column=3, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="m/z tol. for MS2 matches in mDa:").grid(row=0, column=4, sticky="w", padx=(10, 0))
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_ms2_mda, width=10).grid(row=0, column=5, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="metric for scoring MS2 matches").grid(row=0, column=6, sticky="w",
                                                                                 padx=(10, 0))
        self.metric = ttk.Combobox(parameter_frame, values=['cosine_fwd', 'cosine_bwd', 'cosine_sim', 'modified_cosine_greedy'], width=10)
        self.metric.set('modified_cosine_greedy')
        self.metric.grid(row=0, column=7, padx=5,sticky="w")

        ttk.Checkbutton(parameter_frame, text="require_ms2", variable=self.require_ms2).grid(row=0, column=8, padx=(10, 0), sticky="w")
        ttk.Button(parameter_frame, text="Search", command= self.run_find_matches).grid(row=0, column=9, padx=(10, 0), sticky="e")

        #Database status
        self.database_status = ttk.Label(self.root, text="Please Load Database to begin.")
        self.database_status.pack(anchor="w", padx=15, pady=2)

        #Treeview table form
        treeview_frame = ttk.Frame(self.root, padding=10)
        treeview_frame.pack(fill="both", expand=True)

        # this defines the order
        self.tree_columns = ("feature_id_meas",  "name", "formula", "ms2_score", 'n_hits_ms2',
                             "dmz_mda", "dmz_ppm", 'source_library')
        self.tree = ttk.Treeview(treeview_frame, columns=self.tree_columns, show="headings", selectmode="browse")

        #Setting table column headers
        self.tree.heading("feature_id_meas", text="Measurement feature ID")
        self.tree.heading("name", text="Compound name")
        self.tree.heading("formula", text="Formula")
        self.tree.heading("ms2_score", text="MS2 score")
        self.tree.heading("n_hits_ms2", text="#MS2 peak matches")
        self.tree.heading("dmz_mda", text="Delta m/z (mDa)")
        self.tree.heading("dmz_ppm", text="Delta m/z (ppm)")
        self.tree.heading("source_library", text="Source")

        #Setting column widths
        self.tree.column("feature_id_meas", width=60, anchor="center")
        self.tree.column("name", width=250, anchor="w")
        self.tree.column("formula", width=100, anchor="w")
        self.tree.column("ms2_score", width=60, anchor="e")
        self.tree.column("n_hits_ms2", width=60, anchor="center")
        self.tree.column("dmz_mda", width=60, anchor="e")
        self.tree.column("dmz_ppm", width=60, anchor="e")
        self.tree.column("source_library", width=250, anchor="w")

        scroll_bar = ttk.Scrollbar(treeview_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=scroll_bar.set)

        self.tree.pack(side = "left", fill="both", expand=True)
        scroll_bar.pack(side = "right", fill="y")

        button_frame = ttk.Frame(self.root, padding=10)
        button_frame.pack(fill="x", padx=10, pady=5)

        # Plot Area
        self.plot_frame = tk.LabelFrame(self.root, text="Match Overview")
        self.plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        ttk.Button(button_frame, text="Plot Match", command=self.plot_match).pack(side="left", padx=5)
        ttk.Button(button_frame, text = "Annotate feature", command=self.add_annotation_to_meas).pack(side="right", padx=5)

    #Loading Measured Database
    def load_meas_db(self):
        path = filedialog.askopenfilename(filetypes=[("SQLite Files", "*.sqlite"), ("All files", "*.*")])
        if path:
            self.meas_path.set(path)
            try:
                self.meas = FeatureManagerDB(path)
                self.database_status.config(text=f"Measured Database loaded: {path}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load Measured Database:\n {e}")

    #Loading Library Database
    def load_lib_db(self):
        path = filedialog.askopenfilename(filetypes=[("SQLite Files", "*.sqlite"), ("All files", "*.*")])
        if path:
            self.lib_path.set(path)
            try:
                self.lib = Library(path)
                self.database_status.config(text=f"Library Database loaded: {path}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to Load Library Database: \n{e}")

    @functools.cached_property
    def ms2_spectra(self):
        return self.meas.get_ms_spectra(self.meas.feature_ids, level=2)

    def run_find_matches(self):
        if not self.meas_path.get() or not self.lib_path.get():
            messagebox.showerror("Error", "Please select both databases")
            return

        require_ms2: bool = self.require_ms2.get()

        max_dmz_ppm: float = float(self.max_dmz_ppm.get()) if self.max_dmz_ppm.get().strip() else None
        max_dmz_da: float = float(self.max_dmz_mda.get()) * 1e-3 if self.max_dmz_mda.get().strip() else None
        max_ms2_dmz_da: float = float(self.max_dmz_ms2_mda.get()) * 1e-3 if self.max_dmz_ms2_mda.get().strip() else None
        metric = self.metric.get()

        if (max_dmz_da is None and max_dmz_ppm is None) or (max_dmz_ppm is not None and max_dmz_da is not None):
            messagebox.showerror("Error", "Please provide either max_dmz_ppm or max_dmz_mda (but not both).")
            return

        for row in self.tree.get_children():
            self.tree.delete(row)

        try:
            self.current_matches: dict[int, list[dict]] = self.lib.find_matches(
                mzs=self.meas.mzs,
                max_dmz_ppm=max_dmz_ppm,
                max_dmz_da=max_dmz_da,
                max_ms2_dmz_da=max_ms2_dmz_da,
                metric=metric,
                return_nhits_ms2=True,
                ms2_spectra=self.ms2_spectra,
                require_ms2=require_ms2)
        except Exception as e:
            self.current_matches = None
            messagebox.showerror("Error", f"Search failed:\n{e}")
            self.database_status.config(text=f"No matches found.")
            return

        # turn results into dataframe
        series = []
        for f_id_meas, matches in self.current_matches.items():
            for match in matches:
                series.append(pd.Series(name=f_id_meas, data=match))

        self.current_matches_table: pd.DataFrame = (
            pd.concat(series, axis=1).T
            .reset_index(drop=False, names='feature_id_meas')
            .rename(columns={'feature_id': 'feature_id_lib'})
            .astype({
                'feature_id_meas': int,
                'feature_id_lib': int,
                'ms2_score': float,
                'name': str,
                'formula': str,
                'dmz_mda': float,
                'dmz_ppm': float,
                'source_library': str,
                'n_hits_ms2': int}
            )
            .sort_values(by=['feature_id_meas', 'ms2_score'])
        )
        self.current_matches_table.loc[:, 'ms2_score'] = self.current_matches_table['ms2_score'].round(3)
        self.current_matches_table.loc[:, 'dmz_mda'] = self.current_matches_table['dmz_mda'].round(1).abs()
        self.current_matches_table.loc[:, 'dmz_ppm'] = self.current_matches_table['dmz_ppm'].round(1).abs()

        for _, row in self.current_matches_table.iterrows():
            self.tree.insert('', 'end', values=[row[col_name] for col_name in self.tree_columns])

    def plot_match(self):
        selection = self.tree.selection()
        if not selection:
            messagebox.showerror("Error", "Please select a match row first")
            return

        values = self.tree.item(selection, 'values')

        f_id_meas, f_id_lib = int(values[0]), int(values[1])

        match_result = [m for m in self.current_matches[f_id_meas] if m['feature_id'] == f_id_lib][0]

        max_dmz_ppm: float = float(self.max_dmz_ppm.get()) if self.max_dmz_ppm.get().strip() else None
        max_dmz_da: float = float(self.max_dmz_mda.get()) * 1e-3 if self.max_dmz_mda.get().strip() else None
        if max_dmz_da is not None:
            mz_tol_da = max_dmz_da
        else:  # use m/z to convert ppm to da
            mz_tol_da = max_dmz_ppm / 1e6 * self.meas.mzs[f_id_meas] * 1e3

        # Remove the previous canvas
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        if self.toolbar is not None:
            self.toolbar.destroy()
            self.toolbar = None

        fig = Figure(figsize=(15, 15))
        self.lib.plot_match(
            f_id_lib=f_id_lib,
            f_id_meas=f_id_meas,
            meas=self.meas,
            mz_tol_da=mz_tol_da,
            annotation_relative_cutoff=0.3,  # TODO: turn this into slider
            match_result=match_result,
            fig=fig
        )
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.toolbar.pack()
        self.canvas.draw()
        # plt.close(fig)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def add_annotation_to_meas(self):
        selection = self.tree.selection()
        if not selection:
            messagebox.showerror("Error", "Please select a match row first")
            return

        values = self.tree.item(selection, 'values')
        f_id_mas, f_id_lib = int(values[0]), int(values[1])

        confirm = messagebox.askyesno("Confirm", f"Assign {f_id_lib} to Feature {f_id_mas}?")

        if not confirm:
            return
        try:
            add_annotation_to_measured_db(f_id_mas, f_id_lib, self.current_matches_table, self.lib, self.meas)
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = AnnotationGUI(root)
    root.mainloop()
