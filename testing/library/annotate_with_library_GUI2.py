import sqlite3
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import messagebox, filedialog, ttk
import matplotlib.pyplot as plt
from msIO.feature_managers.db import Library, FeatureManagerDB


class AnnotationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Annotation GUI")
        self.root.geometry("950x700")
        self.meas_path = tk.StringVar()
        self.lib_path = tk.StringVar()
        self.f_id = tk.IntVar(value=1)
        self.max_dmz_ppm = tk.StringVar(value="10.0")
        self.max_dmz_da = tk.StringVar(value="")
        self.require_ms2 = tk.BooleanVar(value=False)

        self.meas = None
        self.lib = None
        self.current_matches = None

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

        ttk.Label(parameter_frame, text="Feature ID (f_id):").grid(row=0, column=0, sticky="w")
        ttk.Entry(parameter_frame, textvariable=self.f_id, width=10).grid(row=0, column=1, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="max_dmz_ppm:").grid(row=0, column=2, sticky="w", padx=(10,0))
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_ppm, width=10).grid(row=0, column=3, padx=5, sticky="w")

        ttk.Label(parameter_frame, text="max_dmz_da:").grid(row=0, column=4, sticky="w", padx=(10,0))
        ttk.Entry(parameter_frame, textvariable=self.max_dmz_da, width=10).grid(row=0, column=5, padx=5, sticky="w")

        ttk.Checkbutton(parameter_frame, text="require_ms2", variable=self.require_ms2).grid(row=1, column=0, columnspan=2, pady=5, sticky="w")
        ttk.Button(parameter_frame, text="Search", command= self.run_find_matches).grid(row=1, column=5, sticky="e")

        #Database status
        self.database_status = ttk.Label(self.root, text="Please Load Database to begin.")
        self.database_status.pack(anchor="w", padx=15, pady=2)

        #Treeview table form. Got to know about it from Coders Legacy website
        treeview_frame = ttk.Frame(self.root, padding=10)
        treeview_frame.pack(fill="both", expand=True)

        columns = ("id", "mz","rt","annotation","delta_mz","error_ppm")
        self.tree = ttk.Treeview(treeview_frame, columns=columns, show="headings", selectmode="browse")

        #Setting table column headers
        self.tree.heading("id", text="ID")
        self.tree.heading("mz", text="m/z")
        self.tree.heading("rt", text="RT")
        self.tree.heading("annotation", text="Annotation")
        self.tree.heading("delta_mz", text="Delta m/Z")
        self.tree.heading("error_ppm", text="Error PPM")

        #Setting column widths
        self.tree.column("id", width=60, anchor="center")
        self.tree.column("mz", width=110, anchor="e")
        self.tree.column("rt", width=80, anchor="center")
        self.tree.column("annotation", width=250, anchor="w")
        self.tree.column("delta_mz", width=100, anchor="e")
        self.tree.column("error_ppm", width=100, anchor="e")

        scroll_bar = ttk.Scrollbar(treeview_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=scroll_bar.set)

        self.tree.pack(side = "left", fill="both", expand=True)
        scroll_bar.pack(side = "right", fill="y")

        button_frame = ttk.Frame(self.root, padding=10)
        button_frame.pack(fill="x", padx=10, pady=5)

        ttk.Button(button_frame, text="Plot Match", command=self.plot_match).pack(side="left", padx=5)
        ttk.Button(button_frame, text = "Annotate feature", command=self.annotate_feature).pack(side="right", padx=5)

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

    def run_find_matches(self):
        if not self.meas_path.get() or not self.lib_path.get():
            messagebox.showerror("Error", "Please select both databases")
            return

        f_id = self.f_id.get()
        require_ms2 = self.require_ms2.get()

        max_dmz_ppm = float(self.max_dmz_ppm.get()) if self.max_dmz_ppm.get().strip() else None
        max_dmz_da = float(self.max_dmz_da.get()) if self.max_dmz_da.get().strip() else None

        if (max_dmz_da is None and max_dmz_ppm is None) or (max_dmz_ppm is not None and max_dmz_da is not None):
            messagebox.showerror("Error", "Please provide either max_dmz_ppm or max_dmz_da (but not both).")
            return

        for row in self.tree.get_children():
            self.tree.delete(row)

        matches = None

        try:
            matches = self.lib.find_matches(f_id, self.meas, max_dmz_ppm, max_dmz_da, require_ms2 = require_ms2)

        except TypeError as e:
            if "has no len()" in str(e):
                connection_meas = sqlite3.connect(self.meas_path.get())
                row = connection_meas.execute("SELECT mz FROM peak WHERE id = ?;", (f_id,)).fetchone()
                connection_meas.close()

                if row and row[0] is not None:
                    m_mz = float(row[0])
                    #Calculating lower and upper m/z bounds
                    delta = m_mz * (max_dmz_ppm/1e6) if max_dmz_ppm else max_dmz_da
                    connection_lib = sqlite3.connect(self.lib_path.get())
                    matches = pd.read_sql_query("SELECT id, mz, rt, annotation FROM peak WHERE mz BETWEEN ? AND ?;", connection_lib, params=(m_mz - delta, m_mz + delta))
                    connection_lib.close()
                else:
                    return
            else:
                messagebox.showerror("Error", str(e))
                return

        except Exception as e:
            messagebox.showerror("Error", f"Search failed:\n{e}")

        self.current_matches = matches
        connection_meas = sqlite3.connect(self.meas_path.get())
        peak_info = connection_meas.execute("SELECT mz, rt FROM peak WHERE id = ?;", (f_id,)).fetchone()
        connection_meas.close()

        if peak_info and peak_info[0] is not None:
            m_mz = float(peak_info[0])
            m_rt = f"{peak_info[1]:.2f}" if peak_info[1] is not None else None
            self.database_status.config(text=f"Feature {f_id} | m/z: {m_mz:.5f} | RT: ({m_rt})")

        else:
            m_mz = None

        if matches is not None and not matches.empty:
            for _, r in matches.iterrows():
                library_id = int(r['id']) if 'id' in r else None
                library_mz = float(r['mz']) if 'mz' in r else 0.0
                library_rt = f"{r['rt']:.2f}" if ('rt' in r and pd.notnull(r['rt'])) else None
                annotation = str(r['annotation']) if ('annotation' in r and pd.notnull(r['annotation'])) else None

                if m_mz:
                    dmz = library_mz - m_mz
                    ppm = (dmz / m_mz) * 1e6
                    delta_str, ppm_str = f"{dmz:+.5f}", f"{ppm:+.2f}"
                else:
                    delta_str, ppm_str = None, None

                self.tree.insert('', 'end', values=(library_id, f"{library_mz:.5f}", library_rt, annotation, delta_str, ppm_str))

        else:
            self.database_status.config(text=f"Feature {f_id}: No matches found.")

    def plot_match(self):
        selection = self.tree.selection()
        if not selection:
            messagebox.showerror("Error", "Please select a match row first")
            return

        values = self.tree.item(selection, 'values')
        library_id, library_mz_str, _, annotation = values[0], values[1], values[2], values[3]
        f_id = self.f_id.get()

        connection_meas = sqlite3.connect(self.meas_path.get())
        meas_peak = pd.read_sql_query("SELECT mz, intensity, annotation FROM peak WHERE id = ?;", connection_meas, params=(f_id,))
        connection_meas.close()

        connection_lib = sqlite3.connect(self.lib_path.get())
        lib_peak = pd.read_sql_query("SELECT mz, intensity, annotation FROM peak WHERE id = ?;", connection_lib, params=(library_id,))
        connection_lib.close()

        if meas_peak.empty:
            messagebox.showerror("Error", f"No peaks found for {f_id}")
            return

        m_mz = meas_peak['mz'].values
        m_intensity = (meas_peak['intensity'].values/meas_peak['intensity'].max())*1000 if ('intensity' in meas_peak and meas_peak['intensity'].max()>0) else np.full_like(m_mz, 1000)

        if not lib_peak.empty:
            library_mz = lib_peak['mz'].values
            library_intensity = -(lib_peak['intensity'].values/lib_peak['intensity'].max())*1000 if ('intensity' in lib_peak and lib_peak['intensity'].max()>0) else np.full_like(library_mz, -1000)
        else:
            library_mz = [float(library_mz_str)]
            library_intensity = [-1000]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex = False)

        #Top Plot
        ax1.vlines(m_mz, 0, m_intensity, color='red', linestyle='dashed', lw=1.2)
        ax1.vlines(library_mz, 0, library_intensity, color='blue', linestyle='dashed', lw=1.2)
        ax1.axhline(0, color='black', linestyle='dashed', lw=1)

        for mz, intens, label in zip(m_mz, m_intensity, meas_peak.get('annotation', ['']*len(m_mz))):
            text_str = f"{mz:.4f} {label}".strip()
            ax1.text(mz, intens + 50, text_str, rotation = 90, va='bottom', ha = 'center', fontsize = 10, fontweight = 'bold')

        for mz, intens in zip(library_mz, library_intensity):
            text_str = f"{mz:.4f} {annotation}".strip()
            ax1.text(mz, intens - 50, text_str, rotation = 90, va='top', ha = 'center', fontsize = 10)

        ax1.set_ylabel("Intensity")
        ax1.set_xlabel("m/z")
        ax1.set_ylim(-1150,1150)

        #Bottom Plot
        target_mz = m_mz[0]
        ax2.vlines(m_mz, 0, m_intensity, color='red', linestyle='dashed', lw=1.2)
        ax2.vlines(library_mz, 0, library_intensity, color='blue', linestyle='dashed', lw=1)
        ax2.axhline(0, color='red', linestyle='dashed', lw=1)

        #Setting x-lim to 4 Da above and below to avoid any possible text overlapping.
        ax2.set_xlim(target_mz - 4, target_mz+4)
        ax2.set_ylim(-1150,1150)
        ax2.set_ylabel("Fraction")
        ax2.set_xlabel("mass in Da")

        for mz, intens in zip(m_mz, m_intensity):
            if target_mz - 4 <= mz <= target_mz + 4:
                ax2.text(mz, intens + 50, f"{mz:.4f}", rotation = 90, va='bottom', ha = 'center', fontsize = 10, fontweight = 'bold')

        for mz, intens in zip(library_mz, library_intensity):
            if target_mz - 4 <= mz <= target_mz + 4:
                ax2.text(mz, intens - 50, f"{mz:.4f}", rotation = 90, va='bottom', ha = 'center', fontsize = 10, fontweight = 'bold')

        plt.tight_layout()
        plt.show()

    def annotate_feature(self):
        selection = self.tree.selection()
        if not selection:
            messagebox.showerror("Error", "Please select a match row first")
            return

        values = self.tree.item(selection, 'values')
        annotation = values[3]
        f_id = self.f_id.get()

        confirm = messagebox.askyesno("Confirm", f"Assign {annotation} to Feature {f_id}")

        if not confirm:
            return

        try:
            connection = sqlite3.connect(self.meas_path.get())
            connection.execute("UPDATE peak SET annotation = ? WHERE id = ?;", (annotation, f_id))
            connection.commit()
            connection.close()
            messagebox.showinfo("Success", f"Annotation Feature {f_id}.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = AnnotationGUI(root)
    root.mainloop()

