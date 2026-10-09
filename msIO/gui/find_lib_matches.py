import numpy as np
import pandas as pd
import threading
import tkinter as tk
from tkinter import messagebox, filedialog, ttk

from LipidCalculator import CompoundDict
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from msIO.feature_managers.db import Library, FeatureManagerDB, lib_match_result_to_table


def add_annotation_to_measured_db(
        f_id_meas: int,
        f_id_lib: int,
        results: pd.DataFrame,
        lib: Library,
        meas: FeatureManagerDB
) -> None:
    """Add a matched library entry to the measured database."""
    mask_match: pd.Series = (results.feature_id_meas == f_id_meas) & (results.feature_id_lib == f_id_lib)
    if not mask_match.any():
        ms2_score = float('nan')
    else:
        ms2_score: float = results.loc[
            (results.feature_id_meas == f_id_meas) & (results.feature_id_lib == f_id_lib), 'ms2_score'
        ].squeeze()
    meas.add_annotations_from_library(lib, f_id_meas, f_id_lib, ms2_score)


class AnnotationGUI:
    def __init__(self, root: tk.Tk):
        self.tk_root: tk.Tk = root
        self.tk_root.title("Annotation GUI")
        self.tk_root.geometry("1500x750")
        self.tk_meas_path: tk.StringVar = tk.StringVar()
        self.tk_lib_path: tk.StringVar = tk.StringVar()
        self.tk_max_dmz_ppm: tk.DoubleVar = tk.DoubleVar(value=10)
        self.tk_max_dmz_mda: tk.DoubleVar = tk.DoubleVar(value=3)
        self.tk_max_dmz_ms2_mda: tk.DoubleVar = tk.DoubleVar(value=10)
        self.tk_min_ms2_score: tk.DoubleVar = tk.DoubleVar(value=0)
        self.tk_require_ms2 = tk.BooleanVar(value=False)

        self.connected_measurement_db: FeatureManagerDB = None
        self.connected_library_db: Library = None

        self.current_matches_table: pd.DataFrame = pd.DataFrame()
        self.filtered_results: pd.DataFrame = pd.DataFrame()

        self.canvas: FigureCanvasTkAgg = None
        self.toolbar: NavigationToolbar2Tk = None

        self.measurement_load_number: int = 0

        # colors that will be used in the matching table
        self.feature_colors: list[str] = ["#E5FFCC", "#FFCCCC", "#CCFFCC", "#CCFFFF", "#CCCCFF", "#FFE5CC"]
        self.sort_reverse: dict = {}
        self._create_gui()

    def _create_matching_parameters(self):
        # Parameter Frame
        tk_parameter_frame = ttk.LabelFrame(self.tk_root, text="Matching Parameters", padding=10)
        tk_parameter_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(tk_parameter_frame, text="m/z tol. in ppm:").grid(row=0, column=0, sticky="w", padx=(10, 0))
        ttk.Entry(tk_parameter_frame, textvariable=self.tk_max_dmz_ppm, width=10).grid(row=0, column=1, padx=5,
                                                                                       sticky="w")

        ttk.Label(tk_parameter_frame, text="AND").grid(row=0, column=2, sticky="w", padx=(10, 0))
        ttk.Label(tk_parameter_frame, text="m/z tol. in mDa:").grid(row=0, column=3, sticky="w", padx=(10, 0))
        ttk.Entry(tk_parameter_frame, textvariable=self.tk_max_dmz_mda, width=10).grid(row=0, column=4, padx=5,
                                                                                       sticky="w")

        ttk.Label(tk_parameter_frame, text="m/z tol. for MS2 matches in mDa:").grid(row=1, column=0, sticky="w",
                                                                                    padx=(10, 0))
        ttk.Entry(tk_parameter_frame, textvariable=self.tk_max_dmz_ms2_mda, width=10).grid(row=1, column=1, padx=5,
                                                                                           sticky="w")

        ttk.Label(tk_parameter_frame, text="metric for scoring MS2 matches").grid(row=1, column=2, sticky="w",
                                                                                  padx=(10, 0))
        self.tk_ms2_score_metric: ttk.Combobox = ttk.Combobox(
            tk_parameter_frame,
            state="readonly",
            values=['cosine_fwd', 'cosine_bwd', 'cosine_sim', 'modified_cosine_greedy'],
            width=15
        )
        self.tk_ms2_score_metric.set('modified_cosine_greedy')
        self.tk_ms2_score_metric.grid(row=1, column=3, padx=5, sticky="w")

        ttk.Label(tk_parameter_frame, text="min. MS2 score").grid(row=1, column=4, sticky="w", padx=(10, 0))
        ttk.Entry(tk_parameter_frame, textvariable=self.tk_min_ms2_score, width=10).grid(row=1, column=5, sticky="w",
                                                                                         padx=(10, 0))

        ttk.Checkbutton(tk_parameter_frame, text="require_ms2", variable=self.tk_require_ms2).grid(row=1, column=6,
                                                                                                   padx=(10, 0),
                                                                                                   sticky="w")
        self.tk_initialize_matching = ttk.Button(tk_parameter_frame, text="Find Matches",
                                                 command=self.on_click_find_matches)
        self.tk_initialize_matching.grid(row=6, column=0, padx=(10, 0), sticky="e")

    def _create_tree_view(self):
        # Result table
        tk_treeview_frame: ttk.Frame = ttk.Frame(self.tk_root, padding=10)
        tk_treeview_frame.pack(fill="both", expand=True)

        # this defines the order
        self.tree_columns: list[str] = [
            "feature_id_meas", "feature_id_lib", "RT", "name", "formula", "ms2_score", 'n_hits_ms2',
            "dmz_mda", "dmz_ppm", 'source_library'
        ]
        self.tk_tree_view: ttk.Treeview = ttk.Treeview(
            tk_treeview_frame,
            columns=self.tree_columns,
            show="headings",
            selectmode="browse"
        )

        # Setting table column headers
        self.tk_tree_view.heading("feature_id_meas", text="Measurement feature ID")
        self.tk_tree_view.heading("feature_id_lib", text="Library feature ID")
        self.tk_tree_view.heading("RT", text="RT (min)")
        self.tk_tree_view.heading("name", text="Compound name")
        self.tk_tree_view.heading("formula", text="Formula")
        self.tk_tree_view.heading("ms2_score", text="MS2 score")
        self.tk_tree_view.heading("n_hits_ms2", text="#MS2 peak matches")
        self.tk_tree_view.heading("dmz_mda", text="Delta m/z (mDa)")
        self.tk_tree_view.heading("dmz_ppm", text="Delta m/z (ppm)")
        self.tk_tree_view.heading("source_library", text="Source")

        # Setting column widths
        self.tk_tree_view.column("feature_id_meas", width=40, anchor="center")
        self.tk_tree_view.column("feature_id_lib", width=40, anchor="center")
        self.tk_tree_view.column("RT", width=40, anchor="center")
        self.tk_tree_view.column("name", width=120, anchor="w")
        self.tk_tree_view.column("formula", width=80, anchor="w")
        self.tk_tree_view.column("ms2_score", width=40, anchor="e")
        self.tk_tree_view.column("n_hits_ms2", width=40, anchor="e")
        self.tk_tree_view.column("dmz_mda", width=40, anchor="e")
        self.tk_tree_view.column("dmz_ppm", width=40, anchor="e")
        self.tk_tree_view.column("source_library", width=250, anchor="w")

        y_scroll: ttk.Scrollbar = ttk.Scrollbar(tk_treeview_frame, orient="vertical", command=self.tk_tree_view.yview)
        x_scroll: ttk.Scrollbar = ttk.Scrollbar(tk_treeview_frame, orient="horizontal", command=self.tk_tree_view.xview)
        self.tk_tree_view.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.tk_tree_view.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")

        tk_treeview_frame.rowconfigure(0, weight=1)
        tk_treeview_frame.columnconfigure(0, weight=1)

    def _create_filters(self):
        # Filtering
        tk_filter_frame = ttk.LabelFrame(self.tk_root, text="Filters", padding=8)
        tk_filter_frame.pack(fill="x", padx=10, pady=5)

        row1 = ttk.Frame(tk_filter_frame)
        row1.pack(fill="x", anchor="w")

        row2 = ttk.Frame(tk_filter_frame)
        row2.pack(fill="x", anchor="w")

        row3 = ttk.Frame(tk_filter_frame)
        row3.pack(fill="x", anchor="w")

        ttk.Label(row1, text="meas. feature ID:").grid(row=0, column=0, pady=3, padx=(5, 4), sticky="w")
        self.tk_filter_min_f_id_meas = ttk.Entry(row1, width=10)
        self.tk_filter_min_f_id_meas.grid(row=0, column=1, padx=(0, 0), pady=3)
        ttk.Label(row1, text="-").grid(row=0, column=2, pady=3, padx=(0, 0), sticky="w")
        self.tk_filter_max_f_id_meas = ttk.Entry(row1, width=10)
        self.tk_filter_max_f_id_meas.grid(row=0, column=3, padx=(0, 12), pady=3)

        ttk.Label(row1, text="Compound name substrings:").grid(row=0, column=4, pady=3, padx=(5, 4),
                                                               sticky="w")
        self.tk_filter_name = ttk.Entry(row1, width=20)
        self.tk_filter_name.grid(row=0, column=5, padx=(0, 12), pady=3)

        ttk.Label(row1, text="RT (min)").grid(row=0, column=6, pady=3, padx=(5, 4), sticky="w")
        self.tk_filter_min_rt = ttk.Entry(row1, width=10)
        self.tk_filter_min_rt.grid(row=0, column=7, padx=(0, 0), pady=3)
        ttk.Label(row1, text="-").grid(row=0, column=8, pady=3, padx=(0, 0), sticky="w")
        self.tk_filter_max_rt = ttk.Entry(row1, width=10)
        self.tk_filter_max_rt.grid(row=0, column=9, padx=(0, 0), pady=3)

        ttk.Label(row1, text="Formula must contain").grid(row=0, column=10, pady=3, padx=(5, 4), sticky="w")
        self.tk_filter_subformula = ttk.Entry(row1, width=10)
        self.tk_filter_subformula.grid(row=0, column=11, padx=(0, 0), pady=3)

        ttk.Label(row2, text="min. MS2 score").grid(row=0, column=0, pady=3, padx=(5, 4), sticky="w")
        self.tk_filter_ms2_score = ttk.Entry(row2, width=10)
        self.tk_filter_ms2_score.grid(row=0, column=1, padx=(0, 0), pady=3)

        ttk.Label(row2, text="min. #MS2 peak matches").grid(row=0, column=2, pady=3, padx=(5, 4), sticky="w")
        self.tk_filter_n_ms2_hits = ttk.Entry(row2, width=10)
        self.tk_filter_n_ms2_hits.grid(row=0, column=3, padx=(0, 0), pady=3)

        ttk.Label(row2, text="m/z (Da)").grid(row=0, column=4, pady=3, padx=(5, 4), sticky="w")
        self.tk_filter_min_mz = ttk.Entry(row2, width=10)
        self.tk_filter_min_mz.grid(row=0, column=5, padx=(0, 0), pady=3)
        ttk.Label(row2, text="-").grid(row=0, column=6, pady=3, padx=(0, 0), sticky="w")
        self.tk_filter_max_mz = ttk.Entry(row2, width=10)
        self.tk_filter_max_mz.grid(row=0, column=7, padx=(0, 0), pady=3)

        ttk.Button(row3, text="Update filters", command=self.on_click_apply_filters).grid(
            row=0, column=0, padx=5, pady=3
        )
        ttk.Button(row3, text="Remove Filters", command=self.clear_filters).grid(
            row=0, column=1, padx=5, pady=3
        )

    def _create_gui(self):
        # Selecting database
        tk_database_frame: ttk.LabelFrame = ttk.LabelFrame(self.tk_root, text="Database selection", padding=10)
        tk_database_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(tk_database_frame, text="Measured Database").grid(row=0, column=0, sticky="w")
        # connect to meas_path variable
        ttk.Entry(tk_database_frame, textvariable=self.tk_meas_path, width=55).grid(row=0, column=1, padx=5)
        # button will
        ttk.Button(tk_database_frame, text="Browse", command=self.load_meas_db).grid(row=0, column=2)

        ttk.Label(tk_database_frame, text="Library Database").grid(row=1, column=0, sticky="w")
        ttk.Entry(tk_database_frame, textvariable=self.tk_lib_path, width=55).grid(row=1, column=1, padx=5)
        ttk.Button(tk_database_frame, text="Browse", command=self.on_browse_load_lib).grid(row=1, column=2)

        self._create_matching_parameters()

        self._create_filters()

        self._create_tree_view()

        # Plotting and Annotations button
        tk_button_frame = ttk.Frame(self.tk_root, padding=10)
        tk_button_frame.pack(fill="x", padx=10, pady=3)

        ttk.Button(tk_button_frame, text="Plot Match", command=self.plot_match).pack(side="left", padx=5)
        ttk.Button(tk_button_frame, text="Annotate Feature", command=self.add_annotation_to_meas).pack(side="right",
                                                                                                       padx=5)

        # Plot Area
        self.tk_plot_frame = tk.LabelFrame(self.tk_root, text="Match Overview")
        self.tk_plot_frame.pack(fill="both", expand=True, padx=10, pady=8)

        # Database status
        self.tk_status_field = ttk.Label(self.tk_root, text="Please Load Both Databases to begin.")
        self.tk_status_field.pack(anchor="sw", padx=15, pady=3)

    def clear_tree_view_match_results(self):
        for item in self.tk_tree_view.get_children():
            self.tk_tree_view.delete(item)

    def load_meas_db(self):
        path: str = filedialog.askopenfilename(filetypes=[("SQLite Files", "*.sqlite"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.measurement_load_number += 1
            self.connected_measurement_db: FeatureManagerDB = FeatureManagerDB(path)
            self.tk_meas_path.set(path)  # add file from file dialog to text box

            # Clear data belonging to previous other database
            self.ms2_spectra = None
            self.ms2_loading_error = None
            self.ms2_loading = True

            self.current_matches_table = pd.DataFrame(columns=self.tree_columns)
            self.filtered_results = self.current_matches_table.copy()
            self.clear_tree_view_match_results()

            self.tk_status_field.config(text=f"Measured Database loaded: {path}")

            # Retreiving spectra\
            threading.Thread(
                target=self.load_ms2_spectra_worker,
                args=(self.connected_measurement_db, self.measurement_load_number),
                daemon=True
            ).start()
        except Exception as error:
            self.connected_measurement_db = None
            self.ms2_loading = False
            messagebox.showerror("Error", f"Failed to load Measured Database:\n {error}")

    def load_ms2_spectra_worker(self, meas, load_number):
        try:
            spectra = meas.get_ms_spectra(meas.feature_ids, level=2)
            error_text = None
        except Exception as error:
            spectra = None
            error_text = str(error)

        try:
            self.tk_root.after(0, lambda: self.finish_ms2_loading(load_number, spectra, error_text))
        except tk.TclError:
            pass

    def finish_ms2_loading(self, load_number, spectra, error_text):
        if load_number != self.measurement_load_number:
            return
        self.ms2_loading = False

        if error_text is not None:
            self.ms2_spectra = None
            self.ms2_loading_error = error_text

            self.tk_status_field.config(text=f"Measured Database loaded, but MS2 spectra could not be loaded.")
            messagebox.showerror("MS2 loading error", f"Could not load MS2 spectra:\n {error_text}")
            return
        self.ms2_spectra = spectra
        self.ms2_loading_error = None
        self.tk_status_field.config(text=f"Measured Database and MS2 spectra loaded")

    # Loading Library Database
    def on_browse_load_lib(self):
        path = filedialog.askopenfilename(filetypes=[("SQLite Files", "*.sqlite"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.connected_library_db = Library(path)
            self.tk_lib_path.set(path)
            self.tk_status_field.config(text=f"Library Database loaded: {path}")
        except Exception as error:
            self.connected_library_db = None
            messagebox.showerror("Error", f"Failed to Load Library Database: \n{error}")

    def on_click_find_matches(self):
        if self.connected_measurement_db is None or self.connected_library_db is None:
            messagebox.showerror("Error", "Please select both databases")
            return

        require_ms2: bool = self.tk_require_ms2.get()

        if require_ms2:
            if self.ms2_loading:
                messagebox.showinfo("MS2 spectra are still loading. Try again shortly")
                return
            if self.ms2_loading_error:
                messagebox.showerror("MS2 ERROR", "Failed to load MS2 spectra. Please try again later.")
                return
            if self.ms2_spectra is None:
                messagebox.showwarning("MS2 Not Available")
                return

        max_dmz_ppm: float = self.tk_max_dmz_ppm.get()
        max_dmz_da: float = self.tk_max_dmz_mda.get() * 1e-3
        max_ms2_dmz_da: float = self.tk_max_dmz_ms2_mda.get() * 1e-3
        metric: str = self.tk_ms2_score_metric.get()
        min_ms2_score: float = self.tk_min_ms2_score.get()

        if max_dmz_ppm <= 0:
            messagebox.showerror("Error", "Please enter a positive number.")
            return

        if max_dmz_da <= 0:
            messagebox.showerror("Error", "Please enter a positive number.")
            return

        if not 0 <= min_ms2_score <= 1:
            messagebox.showerror("Invalid Input", "Please enter a MS2 score between 0 and 1.")
            return

        if max_ms2_dmz_da <= 0:
            messagebox.showerror("Invalid MS2 Tolerance", "The number should be greater than 0.")
            return

        self.clear_tree_view_match_results()
        self.tk_status_field.config(text=f"Searching the library...")
        self.tk_initialize_matching.config(state="disabled")
        self.tk_root.update_idletasks()

        try:
            self.current_matches: dict[int, list[dict]] = self.connected_library_db.find_matches(
                mzs=self.connected_measurement_db.mzs,
                max_dmz_ppm=max_dmz_ppm,
                max_dmz_da=max_dmz_da,
                max_ms2_dmz_da=max_ms2_dmz_da,
                metric=metric,
                min_ms2_score=min_ms2_score,
                return_nhits_ms2=True,
                ms2_spectra=self.ms2_spectra,
                require_ms2=require_ms2
            )
        except Exception as e:
            self.current_matches = {}
            messagebox.showerror("Error", f"Search failed:\n{e}")
            self.tk_status_field.config(text=f"Failed while looking for matches.")
            return

        # turn results into dataframe
        self.current_matches_table = lib_match_result_to_table(self.current_matches)
        self.current_matches_table.loc[:, 'ms2_score'] = self.current_matches_table['ms2_score'].round(3)
        self.current_matches_table.loc[:, 'dmz_mda'] = self.current_matches_table['dmz_mda'].round(1).abs()
        self.current_matches_table.loc[:, 'dmz_ppm'] = self.current_matches_table['dmz_ppm'].round(1).abs()
        self.current_matches_table.loc[:, 'RT'] = self.current_matches_table.feature_id_meas.apply(
            lambda f_id: self.connected_measurement_db.retention_times_in_seconds[f_id] / 60
        ).round(2)

        try:
            self.on_click_apply_filters(show_empty_message=False)
            count = len(self.current_matches_table)
            self.tk_status_field.config(text=f"Found {count} matches.")

            if count == 0:
                messagebox.showinfo("Search complete", "No matches found.")
        except Exception as error:
            self.current_matches = None
            self.current_matches_table = pd.DataFrame(columns=self.tree_columns)
            self.filtered_results = self.current_matches_table.copy()
            self.clear_tree_view_match_results()

            messagebox.showerror("Error", f"Search failed:\n{error}")
            self.tk_status_field.config(text=f"No matches found.")

        finally:
            self.tk_initialize_matching.config(state="normal")

    # Filter Results
    def on_click_apply_filters(self, show_empty_message=True):
        if self.current_matches_table.empty:
            self.filtered_results = self.current_matches_table.copy()
            self.display_dataframe(self.filtered_results)
            return

        filtered_matches: pd.DataFrame = self.current_matches_table.copy()

        mask = np.ones(filtered_matches.shape[0], dtype=bool)

        if min_f_id_meas := self.tk_filter_min_f_id_meas.get().strip():
            try:
                min_f_id_meas = int(min_f_id_meas)
                mask &= filtered_matches.feature_id_meas >= min_f_id_meas
            except ValueError:
                messagebox.showerror('Value Error', 'Minimum ID must be a whole number.')
                return

        if max_f_id_meas := self.tk_filter_max_f_id_meas.get().strip():
            try:
                max_f_id_meas = int(max_f_id_meas)
                mask &= filtered_matches.feature_id_meas <= max_f_id_meas
            except ValueError:
                messagebox.showerror('Value Error', 'Maximum ID must be a whole number.')
                return

        if name_substrings := self.tk_filter_name.get().strip():
            for substring in name_substrings.split(';'):
                mask &= filtered_matches.loc[:, "name"].astype(str).str.contains(substring.strip(), case=False, na=False, regex=False)

        if min_rt := self.tk_filter_min_rt.get().strip():
            try:
                min_rt = float(min_rt)
                mask &= filtered_matches.RT >= min_rt
            except ValueError:
                messagebox.showerror('Value Error', 'Minimum RT must be a number.')
        if max_rt := self.tk_filter_max_rt.get().strip():
            try:
                max_rt = float(max_rt)
                mask &= filtered_matches.RT <= max_rt
            except ValueError:
                messagebox.showerror('Value Error', 'Maximum RT must be a whole number.')

        if subformula := self.tk_filter_subformula.get().strip():
            try:
                cd = CompoundDict(subformula)
                print(f'checking for {cd}')
                mask_f = []
                for _, row in filtered_matches.iterrows():
                    cd_m = CompoundDict(row.formula)
                    for el, c in cd.composition.items():
                        if cd_m.composition.get(el, 0) < c:
                            mask_f.append(False)
                            break
                    else:
                        mask_f.append(True)
                mask &= np.array(mask_f)
            except Exception as e:
                messagebox.showerror('Invalid Formula', f"Invalid formula: {e}")
                return

        if min_ms2_score := self.tk_filter_ms2_score.get().strip():
            try:
                min_ms2_score = float(min_ms2_score)
                mask &= filtered_matches.ms2_score >= min_ms2_score
            except ValueError:
                messagebox.showerror('Value Error', 'Minimum MS2 score must be a number.')

        if min_ms2_hits := self.tk_filter_n_ms2_hits.get().strip():
            try:
                min_ms2_hits = int(min_ms2_hits)
                mask &= filtered_matches.n_hits_ms2 >= min_ms2_hits
            except ValueError:
                messagebox.showerror('Value Error', 'Minimum MS2 hits must be a whole number.')

        if min_mz := self.tk_filter_min_mz.get().strip():
            try:
                min_mz = float(min_mz)
                mask &= filtered_matches.feature_id_meas.apply(lambda f_id: self.connected_measurement_db.mzs[f_id]) >= min_mz
            except ValueError:
                messagebox.showerror('Value Error', 'Minimum m/z must be a number.')

        if max_mz := self.tk_filter_max_mz.get().strip():
            try:
                max_mz = float(max_mz)
                mask &= filtered_matches.feature_id_meas.apply(lambda f_id: self.connected_measurement_db.mzs[f_id]) <= max_mz
            except ValueError:
                messagebox.showerror('Value Error', 'Maximum m/z must be a number.')

        self.filtered_results = filtered_matches.loc[mask, :].reset_index(drop=True)
        self.display_dataframe(self.filtered_results)
        if show_empty_message and self.filtered_results.empty:
            self.tk_status_field.config(text=f"No matches meet the current filters.")

        elif not self.filtered_results.empty:
            self.tk_status_field.config(
                text=f"Showing {len(self.filtered_results)}" f"of {len(self.current_matches_table)} matches.")

    def clear_filters(self):
        self.filtered_results = self.current_matches_table.copy()
        self.display_dataframe(self.filtered_results)
        self.tk_status_field.config(text=f"Showing all {len(self.current_matches_table)} matches.")

    # def sort_results(self, column):
    #     if self.filtered_results.empty:
    #         return
    #
    #     if column not in self.filtered_results.columns:
    #         return
    #
    #     reverse = self.sort_reverse.get(column, False)
    #     ascending = not reverse
    #
    #     sort_key = (lambda series: series.astype(str).str.casefold()
    #     if column in ("name", "formula", "source_library")
    #     else series
    #                 )
    #
    #     sorted_results = self.filtered_results.sort_values(by=column, ascending=ascending, na_position="last",
    #                                                        key=sort_key).reset_index(drop=True)
    #     self.sort_reverse[column] = not reverse
    #
    #     self.filtered_results = sorted_results
    #     self.display_dataframe(self.filtered_results)

    def display_dataframe(self, dataframe: pd.DataFrame):
        self.clear_tree_view_match_results()

        # Mapping each measured feature ID to a background color
        f_ids = sorted(dataframe['feature_id_meas'].unique())
        n_colors = len(self.feature_colors)
        f_id_meas_to_bg_color: dict[int, str] = {
            f_id: self.feature_colors[c_idx % n_colors] for c_idx, f_id in enumerate(f_ids)
        }

        for _, row in dataframe.iterrows():
            feature_id_meas: int = row.feature_id_meas

            color: str = f_id_meas_to_bg_color[feature_id_meas]
            tag: str = f"feature_{feature_id_meas}"
            self.tk_tree_view.tag_configure(tag, background=color)

            values: list = []
            for column in self.tree_columns:
                value = row.get(column, "")
                if pd.isna(value):
                    value = ""
                values.append(value)
            self.tk_tree_view.insert("", 'end', values=values, tags=(tag,))

    # Plot Selected Match
    def plot_match(self):
        selected_item: tuple[str, ...] = self.tk_tree_view.selection()
        if not selected_item:
            messagebox.showerror("Error", "Please select a match row first")
            return

        if self.current_matches_table.empty:
            messagebox.showerror("Error", "Run a successful search before plotting a match")
            return

        selected_item: str = selected_item[0]
        # Adding some lines here due to repetitive error for Plot Match
        # Getting measured ID
        try:
            # Finding which row was selected
            selected_index: int = self.tk_tree_view.index(selected_item)
            selected_row: pd.Series = self.filtered_results.iloc[selected_index]
            f_id_meas: int = int(selected_row.loc["feature_id_meas"])
            f_id_lib: int = int(selected_row.loc["feature_id_lib"])
        except (ValueError, IndexError, TypeError, KeyError) as error:
            messagebox.showerror("Error",
                                 f"Could not determine the measured and library feature ID from the selected result:\n{error}")
            return

        match_result = [m for m in self.current_matches[f_id_meas] if m['feature_id'] == f_id_lib][0]

        try:
            max_dmz_ppm: float = self.tk_max_dmz_ppm.get()
            max_dmz_da: float = self.tk_max_dmz_mda.get() * 1e-3

        except ValueError:
            messagebox.showerror("Invalid Tolerance", "Check the precursor tolerance values.")
            return

        mz_tol_da = max(max_dmz_da, max_dmz_ppm * 1e-6 * self.connected_measurement_db.mzs[f_id_meas])
        # Remove the previous canvas
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        if self.toolbar is not None:
            self.toolbar.destroy()
            self.toolbar = None

        fig = Figure(figsize=(15, 15))
        self.connected_library_db.plot_match(
            f_id_lib=f_id_lib,
            f_id_meas=f_id_meas,
            meas=self.connected_measurement_db,
            mz_tol_da=mz_tol_da,
            annotation_relative_cutoff=0.3,  # TODO: turn this into slider
            match_result=match_result,
            fig=fig
        )
        self.canvas = FigureCanvasTkAgg(fig, master=self.tk_plot_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.tk_plot_frame)
        self.toolbar.update()
        self.toolbar.pack()
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def add_annotation_to_meas(self):
        selection = self.tk_tree_view.selection()
        if not selection:
            messagebox.showerror("Error", "Please select a match row first")
            return

        if (
                self.connected_measurement_db is None or self.connected_library_db is None) or self.current_matches_table.empty:
            messagebox.showerror("Error", "Run a successful serch before plotting a match")
            return

        values = self.tk_tree_view.item(selection[0], "values")
        f_id_meas, f_id_lib = int(values[0]), int(values[1])

        confirm = messagebox.askyesno("Confirm", f"Assign {f_id_lib} to Feature {f_id_meas}?")

        if not confirm:
            return
        try:
            add_annotation_to_measured_db(f_id_meas, f_id_lib, self.current_matches_table, self.connected_library_db,
                                          self.connected_measurement_db)
            messagebox.showinfo("Annotation", f"Annotation added for measured feature {f_id_meas}")

        except Exception as error:
            messagebox.showerror("Annotation Error", str(error))


if __name__ == "__main__":
    root: tk.Tk = tk.Tk()
    app: AnnotationGUI = AnnotationGUI(root)
    app.tk_root.mainloop()
