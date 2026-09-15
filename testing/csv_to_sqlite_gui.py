"""Written by Avin Jain"""
import os
import tkinter as tk
from tkinter import filedialog, messagebox
from typing import cast

from msIO.feature_managers.combined import ProjectImportManager
from msIO.feature_managers.metaboscape import MetaboscapeImportManager
from msIO.feature_managers.gnps import GnpsImportManager
from msIO.feature_managers.sirius import SiriusImportManager
from msIO import MgfImportManager
from msIO.sql.session import initiate_db


def select_csv():
    path = filedialog.askopenfilename(
        title="Select CSV file",
        filetypes=[
            ("CSV files", "*.csv"),
            ("All files", "*.*")
        ]
    )

    if path:
        csv_entry.delete(0, tk.END)
        csv_entry.insert(0, path)


def select_mgf():
    path = filedialog.askopenfilename(
        title="Select MGF file",
        filetypes=[
            ("MGF files", "*.mgf"),
            ("All files", "*.*")
        ]
    )

    if path:
        mgf_entry.delete(0, tk.END)
        mgf_entry.insert(0, path)


def select_gnps():
    path = filedialog.askdirectory(
        title="Select GNPS folder"
    )

    if path:
        gnps_entry.delete(0, tk.END)
        gnps_entry.insert(0, path)


def select_sirius():
    path = filedialog.askdirectory(
        title="Select Sirius folder"
    )

    if path:
        sirius_entry.delete(0, tk.END)
        sirius_entry.insert(0, path)


def select_sqlite():
    path = filedialog.asksaveasfilename(
        title="Select SQLite database",
        defaultextension=".sqlite",
        filetypes=[
            ("SQLite files", "*.sqlite"),
            ("All files", "*.*")
        ]
    )

    if path:
        sqlite_entry.delete(0, tk.END)
        sqlite_entry.insert(0, path)


def toggle_mz_resolution():
    if accumulate_var.get():
        mz_resolution_entry.config(state="normal")
    else:
        mz_resolution_entry.config(state="disabled")


def convert():
    path_metaboscape_csv = csv_entry.get()
    path_mgf = mgf_entry.get()
    path_gnps = gnps_entry.get()
    path_sirius = sirius_entry.get()
    export_tag = sirius_tag_entry.get()
    db_file = sqlite_entry.get()

    # MetaboScape
    if not path_metaboscape_csv:
        messagebox.showerror(
            "Error",
            "Please select a CSV file."
        )
        return

    if not db_file:
        messagebox.showerror(
            "Error", "Please select a SQLite database."
        )
        return

    try:
        # MetaboScape
        metaboscape_manager = MetaboscapeImportManager(
            path_metaboscape_csv
        )

        # MGF
        mgf_manager = None
        # Only create one if the user selected an MGF file.
        if path_mgf:
            if accumulate_var.get():
                mz_resolution = mz_resolution_entry.get()
                if not mz_resolution:
                    messagebox.showerror("Error", "Please enter an m/z resolution."
                                         )
                    return

                # Convert the text from the GUI into an integer.
                mz_resolution = int(mz_resolution)
                mgf_manager = MgfImportManager(
                    path_mgf,
                    accumulate_spectra=True,
                    mz_resolution=mz_resolution
                )
            else:
                mgf_manager = MgfImportManager(
                    path_mgf,
                    accumulate_spectra=False
                )

        # GNPS
        gnps_manager = None
        if path_gnps:
            gnps_manager = GnpsImportManager(
                path_gnps_folder=path_gnps
            )

        # SIRIUS
        sirius_manager = None
        if path_sirius:
            if export_tag:
                sirius_manager = SiriusImportManager(
                    path_folder_export=path_sirius,
                    export_tag=export_tag
                )
            else:
                sirius_manager = SiriusImportManager(
                    path_folder_export=path_sirius
                )

        # Combine all import managers
        project_import_manager = ProjectImportManager(
            metaboscape_manager=metaboscape_manager,
            mgf_manager=mgf_manager,
            gnps_manager=gnps_manager,
            sirius_manager=sirius_manager
        )

        # Creating SQLite database
        initiate_db(db_file)
        project_import_manager.to_sql(db_file, feature_ids=metaboscape_manager.feature_ids)

        messagebox.showinfo("Success!", "CSV successfully converted to SQLite.")

    except Exception as e:
        messagebox.showerror("Error!", str(e))


# GUI

root = tk.Tk()
root.title("CSV to SQLite Converter")
root.geometry("750x700")

root.columnconfigure(0,
                     weight=1)  # Allowing the main column to expand when the window is resized (Got to know about this from StackOverflow)

# Metaboscape

metaboscape_frame = tk.LabelFrame(root, text="Metaboscape", pady=10, padx=10)
metaboscape_frame.grid(row=0, column=0, padx=20, pady=(20, 10), sticky="ew")

# Expanding the middle coulumn
metaboscape_frame.columnconfigure(1, weight=1)

tk.Label(metaboscape_frame, text="Input CSV").grid(row=0, column=0, padx=10, pady=5, sticky="w")

csv_entry = tk.Entry(metaboscape_frame)
csv_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")

tk.Button(metaboscape_frame, text="Browse...", command=select_csv).grid(row=0, column=2, padx=10, pady=5)

# MGF

mgf_frame = tk.LabelFrame(root, text="MGF", pady=10, padx=10)
mgf_frame.grid(row=1, column=0, padx=20, pady=10, sticky="ew")

# Expanding the middle coulumn
mgf_frame.columnconfigure(1, weight=1)

tk.Label(mgf_frame, text="MGF File").grid(row=0, column=0, padx=10, pady=5, sticky="w")

mgf_entry = tk.Entry(mgf_frame)
mgf_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")

tk.Button(mgf_frame, text="Browse...", command=select_mgf).grid(row=0, column=2, padx=10, pady=5)

accumulate_var = tk.BooleanVar(value=False)

tk.Checkbutton(
    mgf_frame, text="Accumulate spectra", variable=accumulate_var, command=toggle_mz_resolution
).grid(row=1, column=1, padx=10, pady=5, sticky="w")

# MZ resoltuion

tk.Label(mgf_frame, text="m/z Resolution").grid(row=2, column=0, padx=10, pady=5, sticky="w")

mz_resolution_entry = tk.Entry(mgf_frame, state="disabled")

mz_resolution_entry.grid(row=2, column=1, padx=10, pady=5, sticky="ew")

# GNPS

gnps_frame = tk.LabelFrame(root, text="GNPS", pady=10, padx=10)
gnps_frame.grid(row=2, column=0, padx=20, pady=10, sticky="ew")

# Expanding the middle coulumn
gnps_frame.columnconfigure(1, weight=1)

tk.Label(gnps_frame, text="GNPS Folder").grid(row=0, column=0, padx=10, pady=5, sticky="w")

gnps_entry = tk.Entry(gnps_frame)
gnps_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")

tk.Button(gnps_frame, text="Browse...", command=select_gnps).grid(row=0, column=2, padx=10, pady=5)

# SIRIUS

sirius_frame = tk.LabelFrame(root, text="SIRIUS", pady=10, padx=10)
sirius_frame.grid(row=3, column=0, padx=20, pady=(20, 10), sticky="ew")

# Expanding the middle coulumn
sirius_frame.columnconfigure(1, weight=1)

tk.Label(sirius_frame, text="Input Folder").grid(row=0, column=0, padx=10, pady=5, sticky="w")

sirius_entry = tk.Entry(sirius_frame)
sirius_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")

tk.Button(sirius_frame, text="Browse...", command=select_sirius).grid(row=0, column=2, padx=10, pady=5)

tk.Label(sirius_frame, text="Export Tag").grid(row=1, column=0, padx=10, pady=5, sticky="w")

sirius_tag_entry = tk.Entry(sirius_frame, width=20)
sirius_tag_entry.grid(row=1, column=1, padx=10, pady=5, sticky="ew")

# Output

sqlite_frame = tk.LabelFrame(root, text="Output", pady=10, padx=10)
sqlite_frame.grid(row=4, column=0, padx=20, pady=10, sticky="ew")
sqlite_frame.columnconfigure(1, weight=1)

tk.Label(sqlite_frame, text="SQLite File").grid(row=0, column=0, padx=10, pady=5, sticky="w")

sqlite_entry = tk.Entry(sqlite_frame)
sqlite_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")

tk.Button(sqlite_frame, text="Browse...", command=select_sqlite).grid(row=0, column=2, padx=10, pady=5)

# Convert button
tk.Button(root, text="Convert", command=convert, width=15).grid(row=5, column=0, pady=20)

# Starting GUI
root.mainloop()
