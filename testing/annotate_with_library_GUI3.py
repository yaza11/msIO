import functools
import sqlite3
import pandas as pd
import numpy as np
import threading
import tkinter as tk
from tkinter import messagebox, filedialog, ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from rdkit.sping.colors import black, orange

from msIO.feature_managers.db import Library, FeatureManagerDB
from testing.library.annotate_with_library import f_id_meas


def add_annotation_to_measured_db(f_id_meas: int, f_id_lib: int, results: pd.DataFrame, lib: Library, meas: FeatureManagerDB) -> None:
    """Add a matched library entry to the measured database."""
    selected_score  = results.loc[(results["feature_id_meas"] == f_id_meas)] & (results["feature_id_lib"] == f_id_lib), "ms2_score"

    ms2_score: float = (
        selected_score.loc[0]
        if not selected_score.empty
        else float("nan")
    )

    meas.add_annotations_from_library(lib, f_id_meas, f_id_lib, ms2_score)

class AnnotationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Annotation GUI")
        self.root.geometry("1200x850")
        self.meas_path = tk.StringVar()
        self.lib_path = tk.StringVar()
        self.max_dmz_ppm = tk.StringVar(value="10.0")
        self.max_dmz_mda = tk.StringVar(value="")

        self.metric = tk.StringVar(value="modified_cosine_greedy")
        self.min_ms2_score = tk.StringVar(value=0.3)
        self.max_dmz_ms2_da = tk.StringVar(value=0.01)
        self.require_ms2 = tk.BooleanVar(value=False)

        self.meas = None
        self.lib = None
        self.current_matches = None

        self.current_matches_table = pd.DataFrame()
        self.visible_results = pd.DataFrame()

        self.canvas = None

        self.ms2_spectra = None
        self.ms2_loading = False,
        self.ms2_loading_error = None

        self.measurement_load_number = 0

        self.feature_colors = ["red", "green", "blue", "yellow", "white", "orange"]
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

        ttk.Label(parameter_frame, text="Precusrsor m/z tol. in ppm:").grid(row=0, column=0, sticky="w", padx=(10,0))
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_ppm, width=10).grid(row=0, column=1, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="Precursor m/z tol. in mDa:").grid(row=0, column=2, sticky="w", padx=(10,0))
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_mda, width=10).grid(row=0, column=3, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="Metric").grid(row=0, column=4, sticky="w", padx=(10, 0))
        ttk.Entry(parameter_frame, textvariable=self.metric, width=27).grid(row=0, column=5, padx=5)

        ttk.Label(parameter_frame, text="Minimum MS2 score:").grid(row=1, column=0, sticky="w")
        ttk.Entry(parameter_frame, textvariable=self.min_ms2_score, width=10).grid(row=1, column=1, padx=5)

        ttk.Label(parameter_frame, text="MS2 fragment tolerance (Da):").grid(row=1, column=2, sticky="w")
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_ms2_da, width=10).grid(row=1, column=3, padx=5, sticky="w")

        ttk.Checkbutton(parameter_frame, text="require_ms2", variable=self.require_ms2).grid(row=1, column=4, padx=5, sticky="w")
        self.search_button = ttk.Button(parameter_frame, text="Find Matches", command=self.run_find_matches)
        self.search_button.grid(row=1, column=5, padx=5)

        #Database status
        self.database_status = ttk.Label(self.root, text="Please Load Both Database to begin.")
        self.database_status.pack(anchor="w", padx=15, pady=3)

        #Filtering
        filter_frame = ttk.LabelFrame(self.root, text="Filter Results", padding=8)
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

        self.columns = ("feature_id_meas", "feature_id_lib", "name","formula", "ms2_score", "dmz_da", "dmz_ppm", 'source_library', 'n_hits_ms2')
        self.tree = ttk.Treeview(treeview_frame, columns=self.columns, show="headings", selectmode="browse")

        headings = {"feature_id_meas": "Measurement feature ID",
                    "feature_id_lib": "Lib feature ID",
                    "name": "Compound Name",
                    "formula": "Formula",
                    "ms2_score": "MS2 score",
                    "dmz_da": "DMZ da",
                    "dmz_ppm": "DMZ ppm",
                    "source_library": "Source",
                    "n_hits_ms2": "n_hits_ms2"
        }

        for column, label, in headings.items():
            self.tree.heading(
                column,
                text=label,
                command = lambda col = column:
                    self.sort_results(col)
            )

        widths = {
            "feature_id_meas": 145,
            "feature_id_lib": 95,
            "name": 210,
            "formula": 110,
            "ms2_score": 90,
            "dmz_da": 105,
            "dmz_ppm": 105,
            "source_library": 180,
            "n_hits_ms2": 110
        }

        for column, width in widths.items():
            self.tree.column(
                column,
                width=width,
                anchor=("w" if column in ("name", "formula","source_library") else "center")
            )



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
            self.meas = FeatureManagerDB(path)
            self.meas_path.set(path)

            #Clear data belonging to previous other database
            self.ms2_spectra = None
            self.ms2_loading_error = None
            self.ms2_loading = True

            self.current_matches = None
            self.current_matches_table = pd.DataFrame(columns=self.columns)
            self.visible_results = self.current_matches_table.copy()
            self.clear_tree()

            self.database_status.config(text=f"Measured Database loaded: {path}")

            #Retreiving spectra\
            threading.Thread(target=self.load_ms2_spectra_worker, args = (self.meas, load_number), daemon=True).start()


        except Exception as error:
            self.meas = None
            self.ms2_loading = False
            messagebox.showerror("Error", f"Failed to load Measured Database:\n {e}")

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
            messagebox.showerror("Error", f"Failed to Load Library Database: \n{e}")



    def run_find_matches(self):
        if self.meas is None or self.lib is None:
            messagebox.showerror("Error", "Please select both databases")
            return

        # f_id = self.f_id.get()
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
            max_dmz_ppm = (float(self.max_dmz_ppm.get()) if self.max_dmz_ppm.get().strip() else None)
            max_dmz_da = (float(self.max_dmz_mda.get()) * 1e-3 if self.max_dmz_mda.get().strip() else None)
            metric = self.metric.get().strip()
            min_ms2_score = float(self.min_ms2_score.get())
            max_dmz_ms2_da = float(self.max_dmz_ms2_da.get())

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

        if not metric:
            messagebox.showerror("Error", "Please enter a matching metric.")
            return

        if not 0 <= min_ms2_score <= 1:
            messagebox.showerror("Invalid Input", "Please enter a MS2 score between 0 and 1.")
            return

        if max_dmz_ms2_da <= 0:
            messagebox.showerror("Invalid MS2 Tolerance", "The number should be greater than 0.")
            return

        self.clear_tree()
        self.database_status.config(text=f"Searching the library...")
        self.search_button.config(state="disabled")
        self.root.update_idletasks()

        try:
            self.current_matches = self.lib.find_matches(
                mzs=self.meas.mzs,
                max_dmz_ppm=max_dmz_ppm,
                max_dmz_da=max_dmz_da,
                metric = metric,
                ms2_spectra=self.ms2_spectra if require_ms2 else {},
                require_ms2=require_ms2,
                return_nhits_ms2 = True,
                min_ms2_score=min_ms2_score
            )

            #Converting matches into table rows.
            rows = []
            for f_id_meas, matches in self.current_matches.items():
                for match in matches:
                    result_row = dict(match)
                    result_row['feature_id_meas'] = int(f_id_meas)

                    if "feature_id" in result_row:
                        result_row['feature_id'] = (result_row.pop('feature_id'))

                    elif "feature_id_lib" in result_row:
                        result_row['feature_id_lib'] = (result_row.pop('feature_id_lib'))

                    rows.append(result_row)

                if rows:
                    self.current_matches_table = pd.DataFrame(rows)

                    for column in self.columns:
                        if column not in self.current_matches_table.columns:
                            self.current_matches_table[column] = pd.NA

                    for column in ("feature_id_meas", "feature_id_lib", "ms2_score"):
                        self.current_matches_table[column] = pd.to_numeric(self.current_matches_table[column], errors = "coerce")

                    self.current_matches_table = (
                        self.current_matches_table.sort_values(by=["feature_id_meas", "ms2_score"],
                        ascending=[True, False], na_position="last").reset_index(drop=True)
                    )

                else:
                    self.current_matches_table = pd.DataFrame(columns=self.columns)

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
        self.visible_dataframe(self.visible_results)
        self.database_status.config(text=f"Showing all {len(self.current_matches_table)} matches.")

    #Sorting
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

    #Displaying results and color rows
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

            for column in self.columns:
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
        selected_item = self.tree.selection()
        if not selected_item:
            messagebox.showerror("Error", "Please select a match row first")
            return

        if (self.current_matches is None or self.meas is None or self.lib is None):
            messagebox.showerror("Error", "Run a successful serch before plotting a match")
            return

        selected_item = selected_item[0]
        # ("feature_id_meas", "feature_id_lib", "name", "formula", "ms2_score", "dmz_da", "dmz_ppm", 'source_library',
        #  'n_hits_ms2')
        #Adding some lines here due to repetitive error for Plot Match
        #Getting measured ID
        try:
            #Finding which row was selected
            selected_index = self.tree.index(selected_item)
            selected_row = self.visible_results.iloc[selected_index]
            f_id_meas = int(selected_row["feature_id_meas"])
            f_id_lib = int(selected_row["feature_id_lib"])
        except (ValueError, IndexError, TypeError, KeyError):
            messagebox.showerror("Error", "Could not determine the measured and library feature ID from the selected result")
            return

        match_result = next(
            (m for m in self.current_matches.get(f_id_meas, [])
             if int (m.get("feature_id", -1)) == f_id_lib),
            None
        )

        if match_result is None:
            messagebox.showerror("Error", "Could locate the selected match")
            return

        try:
            max_dmz_ppm: float = float(self.max_dmz_ppm.get()) if self.max_dmz_ppm.get().strip() else None
            max_dmz_da: float = float(self.max_dmz_mda.get()) * 1e-3 if self.max_dmz_mda.get().strip() else None

        except ValueError:
            messagebox.showerror("Invalid Tolerance", "Check the precusror tolerance values.")
            return

        if max_dmz_da is not None:
            mz_tol_da = max_dmz_da

        elif max_dmz_ppm is not None:
            mz_tol_da = max_dmz_ppm / 1e6 * self.meas.mzs[f_id_meas] * 1e3

        else:
            messagebox.showerror("Invalid Tolerance", "Enter either ppm or Da values.")

        # Remove the previous canvas
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None

        try:
            fig = plt.figure(figsize=(15, 15))
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
            plt.close(fig)

            self.canvas.draw()
            self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        except Exception as error:
            plt.close("all")
            messagebox.showerror("Plot Error", f"Could not plot the selected match. {error}")

    def add_annotation_to_meas(self):
        selection = self.tree.selection()
        if not selection:
            messagebox.showerror("Error", "Please select a match row first")
            return

        if (self.meas is None or self.lib is None or self.current_matches_table.empty):
            messagebox.showerror("Error", "Run a successful serch before plotting a match")
            return

        values = self.tree.item(selection[0], "values")
        f_id_mas, f_id_lib = int(values[0]), int(values[1])

        confirm = messagebox.askyesno("Confirm", f"Assign {f_id_lib} to Feature {f_id_meas}?")

        if not confirm:
            return
        try:
            add_annotation_to_measured_db(f_id_meas, f_id_lib, self.current_matches_table, self.lib, self.meas)
            messagebox.showinfo("Annotation", f"Annotation added for measured feature {f_id_meas}")

        except Exception as error:
            messagebox.showerror("Annotation Error", str(error))

if __name__ == "__main__":
    root = tk.Tk()
    app = AnnotationGUI(root)
    root.mainloop()

