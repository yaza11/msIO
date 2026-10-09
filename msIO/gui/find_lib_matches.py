import pandas as pd
import threading
import tkinter as tk
from tkinter import messagebox, filedialog, ttk
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
    def __init__(self, root):
        self.root = root
        self.root.title("Annotation GUI")
        self.root.geometry("1500x750")
        self.meas_path = tk.StringVar()
        self.lib_path = tk.StringVar()
        self.max_dmz_ppm = tk.StringVar(value="10.0")
        self.max_dmz_mda = tk.StringVar(value="")
        self.max_dmz_ms2_mda = tk.StringVar(value="10")
        self.min_ms2_score = tk.StringVar(value="0")
        self.require_ms2 = tk.BooleanVar(value=False)

        self.meas: FeatureManagerDB = None
        self.lib = None
        self.current_matches = None

        self.current_matches_table = pd.DataFrame()
        self.visible_results = pd.DataFrame()

        self.canvas = None
        self.toolbar = None

        self.measurement_load_number = 0

        self.feature_colors = ["#FFCCCC", "#E5FFCC", "#CCFFCC", "#CCFFFF", "#CCCCFF", "#FFE5CC"]
        self.sort_reverse = {}
        self.create_gui()

    def create_gui(self):

        #Selecting database
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

        ttk.Label(parameter_frame, text="m/z tol. for MS2 matches in mDa:").grid(row=1, column=0, sticky="w", padx=(10, 0))
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_ms2_mda, width=10).grid(row=1, column=1, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="metric for scoring MS2 matches").grid(row=1, column=2, sticky="w",
                                                                                 padx=(10, 0))
        self.metric = ttk.Combobox(parameter_frame, state="readonly", values=['cosine_fwd', 'cosine_bwd', 'cosine_sim', 'modified_cosine_greedy'], width=10)
        self.metric.set('modified_cosine_greedy')
        self.metric.grid(row=1, column=3, padx=5,sticky="w")

        ttk.Label(parameter_frame, text="min. MS2 score").grid(row=1, column=4, sticky="w", padx=(10, 0))
        ttk.Entry(parameter_frame, textvariable=self.min_ms2_score, width=10).grid(row=1, column=5, sticky="w",
                                                                                 padx=(10, 0))

        ttk.Checkbutton(parameter_frame, text="require_ms2", variable=self.require_ms2).grid(row=1, column=6, padx=(10, 0), sticky="w")
        self.search_button = ttk.Button(parameter_frame, text="Find Matches", command= self.run_find_matches)
        self.search_button.grid(row=6, column=0, padx=(10, 0), sticky="e")

        #Database status
        self.database_status = ttk.Label(self.root, text="Please Load Both Databases to begin.")
        self.database_status.pack(anchor="w", padx=15, pady=3)

        #Filtering
        filter_frame = ttk.LabelFrame(self.root, text="Filters", padding=8)
        filter_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(filter_frame, text="Minimum measured ID:").grid(row=0, column=0, pady=3, padx=(5,4), sticky="w")
        self.min_id_entry = ttk.Entry(filter_frame, width=10)
        self.min_id_entry.grid(row=0, column=1, padx=(0,12), pady=3)

        ttk.Label(filter_frame, text="Maximum measured ID:").grid(row=0, column=2, pady=3, padx=(5,4), sticky = "w")
        self.max_id_entry = ttk.Entry(filter_frame, width=10)
        self.max_id_entry.grid(row=0, column=3, padx=(0,12), pady=3)

        ttk.Label(filter_frame, text="Compound name contains:").grid(row=0, column=4, pady=3, padx=(5,4), sticky="w")
        self.name_filter_entry = ttk.Entry(filter_frame, width=20)
        self.name_filter_entry.grid(row=0, column=5, padx=(0,12), pady=3)

        ttk.Button(filter_frame, text="Apply filters", command = self.apply_filters).grid(row=0, column=6, padx=5, pady=3)
        ttk.Button(filter_frame, text="Clear Filters", command=self.clear_filters).grid(row=0, column=7, padx=5, pady=3)

        #Result table
        treeview_frame = ttk.Frame(self.root, padding=10)
        treeview_frame.pack(fill="both", expand=True)

        # this defines the order
        self.tree_columns = ("feature_id_meas",  "feature_id_lib", "name", "formula", "ms2_score", 'n_hits_ms2',
                             "dmz_mda", "dmz_ppm", 'source_library')
        self.tree = ttk.Treeview(treeview_frame, columns=self.tree_columns, show="headings", selectmode="browse")

        #Setting table column headers
        self.tree.heading("feature_id_meas", text="Measurement feature ID")
        self.tree.heading("feature_id_lib", text="Library feature ID")
        self.tree.heading("name", text="Compound name")
        self.tree.heading("formula", text="Formula")
        self.tree.heading("ms2_score", text="MS2 score")
        self.tree.heading("n_hits_ms2", text="#MS2 peak matches")
        self.tree.heading("dmz_mda", text="Delta m/z (mDa)")
        self.tree.heading("dmz_ppm", text="Delta m/z (ppm)")
        self.tree.heading("source_library", text="Source")

        #Setting column widths
        self.tree.column("feature_id_meas", width=40, anchor="center")
        self.tree.column("feature_id_lib", width=40, anchor="center")
        self.tree.column("name", width=120, anchor="w")
        self.tree.column("formula", width=80, anchor="w")
        self.tree.column("ms2_score", width=40, anchor="e")
        self.tree.column("n_hits_ms2", width=40, anchor="e")
        self.tree.column("dmz_mda", width=40, anchor="e")
        self.tree.column("dmz_ppm", width=40, anchor="e")
        self.tree.column("source_library", width=250, anchor="w")

        y_scroll = ttk.Scrollbar(treeview_frame, orient="vertical", command=self.tree.yview)
        x_scroll = ttk.Scrollbar(treeview_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")

        treeview_frame.rowconfigure(0, weight=1)
        treeview_frame.columnconfigure(0, weight=1)

        #Plotting and Annotations button
        button_frame = ttk.Frame(self.root, padding=10)
        button_frame.pack(fill="x", padx=10, pady=3)

        ttk.Button(button_frame, text="Plot Match", command = self.plot_match).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Annotate Feature", command = self.add_annotation_to_meas).pack(side="right", padx=5)

        # Plot Area
        self.plot_frame = tk.LabelFrame(self.root, text="Match Overview")
        self.plot_frame.pack(fill="both", expand=True, padx=10, pady=8)

    #Loading Measured Database
    def load_meas_db(self):
        path = filedialog.askopenfilename(filetypes=[("SQLite Files", "*.sqlite"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.measurement_load_number +=1
            load_number = self.measurement_load_number
            self.meas: FeatureManagerDB = FeatureManagerDB(path)
            self.meas_path.set(path)

            #Clear data belonging to previous other database
            self.ms2_spectra = None
            self.ms2_loading_error = None
            self.ms2_loading = True

            self.current_matches = None
            self.current_matches_table = pd.DataFrame(columns=self.tree_columns)
            self.visible_results = self.current_matches_table.copy()
            self.clear_tree()

            self.database_status.config(text=f"Measured Database loaded: {path}")

            #Retreiving spectra\
            threading.Thread(target=self.load_ms2_spectra_worker, args = (self.meas, load_number), daemon=True).start()
        except Exception as error:
            self.meas = None
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
            self.root.after(0, lambda: self.finish_ms2_loading(load_number, spectra, error_text))
        except tk.TclError: #Got to know wbout this from Stackoverflow
            pass

    def finish_ms2_loading(self, load_number, spectra, error_text):
        if load_number != self.measurement_load_number:
            return
        self.ms2_loading = False

        if error_text is not None:
            self.ms2_spectra = None
            self.ms2_loading_error = error_text

            self.database_status.config(text=f"Measured Database loaded, but MS2 spectra could not be loaded.")
            messagebox.showerror("MS2 loading error", f"Could not load MS2 spectra:\n {error_text}")
            return
        self.ms2_spectra = spectra
        self.ms2_loading_error = None
        self.database_status.config(text=f"Measured Database and MS2 spectra loaded")

    #Loading Library Database
    def load_lib_db(self):
        path = filedialog.askopenfilename(filetypes=[("SQLite Files", "*.sqlite"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.lib = Library(path)
            self.lib_path.set(path)
            self.database_status.config(text=f"Library Database loaded: {path}")
        except Exception as error:
            self.lib = None
            messagebox.showerror("Error", f"Failed to Load Library Database: \n{error}")

    def run_find_matches(self):
        if self.meas is None or self.lib is None:
            messagebox.showerror("Error", "Please select both databases")
            return

        require_ms2: bool = self.require_ms2.get()

        if require_ms2:
            if self.ms2_loading:
                messagebox.showinfo("MS2 spectra are still loading. Try again shortly")
                return
            if self.ms2_loading_error:
                messagebox.showerror("MS2 ERROR","Failed to load MS2 spectra. Please try again later.")
                return
            if self.ms2_spectra is None:
                messagebox.showwarning("MS2 Not Available")
                return

        try:
            max_dmz_ppm: float = float(self.max_dmz_ppm.get()) if self.max_dmz_ppm.get().strip() else None
            max_dmz_da: float = float(self.max_dmz_mda.get()) * 1e-3 if self.max_dmz_mda.get().strip() else None
            max_ms2_dmz_da: float = float(self.max_dmz_ms2_mda.get()) * 1e-3 if self.max_dmz_ms2_mda.get().strip() else None
            metric = self.metric.get()
            min_ms2_score: float = float(self.min_ms2_score.get()) if self.min_ms2_score.get().strip() else None
        except ValueError:
            messagebox.showerror("Invalid input", "Please enter a numeric value.")
            return

        if (max_dmz_da is None and max_dmz_ppm is None) or (max_dmz_ppm is not None and max_dmz_da is not None):
            messagebox.showerror("Error", "Please provide either max_dmz_ppm or max_dmz_mda (but not both).")
            return

        if max_dmz_ppm is not None and max_dmz_ppm <= 0:
            messagebox.showerror("Error", "Please enter a positive number.")
            return

        if max_dmz_da is not None and max_dmz_da <= 0:
            messagebox.showerror("Error", "Please enter a positive number.")
            return

        if not 0 <= min_ms2_score <= 1:
            messagebox.showerror("Invalid Input", "Please enter a MS2 score between 0 and 1.")
            return

        if max_ms2_dmz_da <= 0:
            messagebox.showerror("Invalid MS2 Tolerance", "The number should be greater than 0.")
            return

        self.clear_tree()
        self.database_status.config(text=f"Searching the library...")
        self.search_button.config(state="disabled")
        self.root.update_idletasks()

        try:
            self.current_matches: dict[int, list[dict]] = self.lib.find_matches(
                mzs=self.meas.mzs,
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
            self.database_status.config(text=f"Failed while looking for matches.")
            return

        # turn results into dataframe
        self.current_matches_table = lib_match_result_to_table(self.current_matches)
        self.current_matches_table.loc[:, 'ms2_score'] = self.current_matches_table['ms2_score'].round(3)
        self.current_matches_table.loc[:, 'dmz_mda'] = self.current_matches_table['dmz_mda'].round(1).abs()
        self.current_matches_table.loc[:, 'dmz_ppm'] = self.current_matches_table['dmz_ppm'].round(1).abs()

        try:
            self.apply_filters(show_empty_message=False)
            count = len(self.current_matches_table)
            self.database_status.config(text=f"Found {count} matches.")

            if count == 0:
                messagebox.showinfo("Search complete", "No matches found.")
        except Exception as error:
            self.current_matches = None
            self.current_matches_table = pd.DataFrame(columns=self.columns)
            self.visible_results = self.current_matches_table.copy()
            self.clear_tree()

            messagebox.showerror("Error", f"Search failed:\n{error}")
            self.database_status.config(text=f"No matches found.")

        finally:
            self.search_button.config(state="normal")

    #Filter Results
    def apply_filters(self, show_empty_message=True):
        if self.current_matches_table.empty:
            self.visible_results = self.current_matches_table.copy()
            self.display_dataframe(self.visible_results)
            return

        filtered = self.current_matches_table.copy()
        try:
            min_text = self.min_id_entry.get().strip()
            max_text = self.max_id_entry.get().strip()

            if min_text:
                filtered = filtered[filtered['feature_id_meas'] >= int(min_text)]

            if max_text:
                filtered = filtered[filtered['feature_id_meas'] <= int(max_text)]

            if (min_text and max_text and int(min_text) > int(max_text)):
                messagebox.showerror("Error", "Minimum ID cannot exceed Maximum ID.")
                return

        except ValueError:
            messagebox.showerror("Error", "Minimum and Maximum ID must be a whole number.")
            return

        name_text = self.name_filter_entry.get().strip()
        if name_text:
            filtered = filtered[filtered["name"].astype(str).str.contains(name_text, case=False, na=False, regex=False)]

        self.visible_results = filtered.reset_index(drop=True)
        self.display_dataframe(self.visible_results)
        if show_empty_message and self.visible_results.empty:
            self.database_status.config(text=f"No matches meet the current filters.")

        elif not self.visible_results.empty:
            self.database_status.config(text=f"Showing {len(self.visible_results)}" f"of {len(self.current_matches_table)} matches.")

    def clear_filters(self):
        self.min_id_entry.delete(0, tk.END)
        self.max_id_entry.delete(0, tk.END)
        self.name_filter_entry.delete(0, tk.END)
        self.visible_results = (self.current_matches_table.copy()).copy()
        self.display_dataframe(self.visible_results)
        self.database_status.config(text=f"Showing all {len(self.current_matches_table)} matches.")

    def sort_results(self, column):
        if self.visible_results.empty:
            return

        if column not in self.visible_results.columns:
            return

        reverse = self.sort_reverse.get(column, False)
        ascending = not reverse

        sort_key = (lambda series: series.astype(str).str.casefold()
        if column in ("name", "formula", "source_library")
        else series
                    )

        sorted_results = self.visible_results.sort_values(by=column, ascending=ascending, na_position="last", key=sort_key).reset_index(drop=True)
        self.sort_reverse[column] = not reverse

        self.visible_results = sorted_results
        self.display_dataframe(self.visible_results)

    def clear_tree(self):
        for item in self.tree.get_children():
            self.tree.delete(item)

    def display_dataframe(self, dataframe):
        self.clear_tree()

        #Mapping each measured feature ID to a background color
        color_by_feature = {}

        for _, row in dataframe.iterrows():
            feature_value = row.get("feature_id_meas")
            if pd.isna(feature_value):
                feature_key = "unknown"
                display_id = ""

            else:
                feature_key = str(int(feature_value))
                display_id = int(feature_value)

            if feature_key not in color_by_feature:
                color_index = (len(color_by_feature) % len(self.feature_colors))
                color_by_feature[feature_key] = (self.feature_colors[color_index])

            tag = f"feature_{feature_key}"
            self.tree.tag_configure(tag, background=color_by_feature[feature_key])

            values = []

            for column in self.tree_columns:
                value = row.get(column, "")
                if pd.isna(value):
                    value = ""

                elif (column == "feature_id_meas" and display_id != ""):
                    value = display_id

                elif column == "feature_id_lib" and value != "":
                    try:
                        value = int(value)
                    except (TypeError, ValueError):
                        pass
                elif column == "ms2_score" and value != "":
                    try:
                        value = f"{float(value):.4f}"
                    except (TypeError, ValueError):
                        pass

                values.append(value)
            self.tree.insert("", 'end', values=values, tag=(tag,))

    #Plot Selected Match
    def plot_match(self):
        selected_item: tuple[str, ...] = self.tree.selection()
        if not selected_item:
            messagebox.showerror("Error", "Please select a match row first")
            return

        if (self.current_matches is None or self.meas is None or self.lib is None):
            messagebox.showerror("Error", "Run a successful serch before plotting a match")
            return

        selected_item: str = selected_item[0]
        # Adding some lines here due to repetitive error for Plot Match
        # Getting measured ID
        try:
            #Finding which row was selected
            selected_index: int = self.tree.index(selected_item)
            selected_row: pd.Series = self.visible_results.iloc[selected_index]
            f_id_meas: int = int(selected_row.loc["feature_id_meas"])
            f_id_lib: int = int(selected_row.loc["feature_id_lib"])
        except (ValueError, IndexError, TypeError, KeyError) as error:
            messagebox.showerror("Error", f"Could not determine the measured and library feature ID from the selected result:\n{error}")
            return

        match_result = [m for m in self.current_matches[f_id_meas] if m['feature_id'] == f_id_lib][0]

        try:
            max_dmz_ppm: float = float(self.max_dmz_ppm.get()) if self.max_dmz_ppm.get().strip() else None
            max_dmz_da: float = float(self.max_dmz_mda.get()) * 1e-3 if self.max_dmz_mda.get().strip() else None

        except ValueError:
            messagebox.showerror("Invalid Tolerance", "Check the precursor tolerance values.")
            return

        if max_dmz_da is not None:
            mz_tol_da = max_dmz_da
        else:  # use m/z to convert ppm to da
            mz_tol_da = max_dmz_ppm / 1e6 * self.meas.mzs[f_id_meas]

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
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def add_annotation_to_meas(self):
        selection = self.tree.selection()
        if not selection:
            messagebox.showerror("Error", "Please select a match row first")
            return

        if (self.meas is None or self.lib is None) or self.current_matches_table.empty:
            messagebox.showerror("Error", "Run a successful serch before plotting a match")
            return

        values = self.tree.item(selection[0], "values")
        f_id_meas, f_id_lib = int(values[0]), int(values[1])

        confirm = messagebox.askyesno("Confirm", f"Assign {f_id_lib} to Feature {f_id_meas}?")

        if not confirm:
            return
        try:
            add_annotation_to_measured_db(f_id_meas, f_id_lib, self.current_matches_table, self.lib, self.meas)
            messagebox.showinfo("Annotation", f"Annotation added for measured feature {f_id_meas}")

        except Exception as error:
            messagebox.showerror("Annotation Error", str(error))

if __name__ == "__main__":
    root: tk.Tk = tk.Tk()
    app: AnnotationGUI = AnnotationGUI(root)
    app.root.mainloop()
