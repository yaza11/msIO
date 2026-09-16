"""Written by Avin Jain"""
from msIO.feature_managers.db import Library
import tkinter as tk
from tkinter import filedialog, messagebox

# #Location of file
# db_file = r"D:/avin college/University of Bremen/Hiwi work/Hiwi with Yannick/msIO/testing/test data/lib/archlipids_high_conf.sqlite"
#
# #Checking tools present inside Library object
# library = Library(db_file)
# print(dir(library))

library = None

#Selecting the SQLite file
def select_database():
    path = filedialog.askopenfilename(
        title="Select SQLite database",
        filetypes=[("sqlite files", "*.sqlite"),
                   ("SQLite databse", "*.db"),
                   ("All files", "*.*")
        ]
    )

    #After user selected a file
    if path:
        db_entry.delete(0, tk.END)
        db_entry.insert(0, path)

#Loading databse
def load_database():

    global library

    db_file=db_entry.get()
    if not db_file:
        messagebox.showerror("Error", "No database selected")
        return

    try:
        library = Library(db_file)
        status_label.config(text="Database loaded successfully")

        messagebox.showinfo("Library loaded", "Library loaded successfully")

    except Exception as e:

        #Show error if something goes wrong
        library = None

        status_label.config(text="Database could not be loaded")

        messagebox.showerror("Error", str(e))

def search_compound():

    #Making sure the database has been loaded
    if library is None:
        messagebox.showerror("Error", "No database selected, please load a SQLite database first")
        return

    compound_name = compound_entry.get().strip()

    if not compound_name:
        messagebox.showerror("Error", "No compound name selected, please select a compound")
        return
    try:
        matches = library.find_by_name(compound_name)

        print("Matches")
        print(matches)

        #Displaying result
        result_text.delete("1.0", tk.END)

        result_text.insert(tk.END, f"Matches:\n{matches}\n")

    except Exception as e:
        messagebox.showerror("Error", str(e))

#Displaying library info
def display_library_info():

    if library is None:
        messagebox.showerror("Error", "No database selected")
        return

    try:
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, "Library information\n")

        #Feaure IDs
        result_text.insert(tk.END, f"Number of features: {len(library.feature_ids)}\n\n")

        #Names
        result_text.insert(tk.END, f"Number of features with names: {len(library.names)}\n\n")

        #First 10 names
        result_text.insert(tk.END, "First available names:\n")

        for idx, name in enumerate(library.names.values()):
            if idx >= 10:
                break
            result_text.insert(tk.END, f"{name}\n")

    except Exception as e:
        messagebox.showerror("Error", str(e))

#Displaying Compund Properties of the searched compound

def show_compound_properties():
    if library is None:
        messagebox.showerror("Error", "No database selected, please load a SQLite database first")
        return

    compound_name = compound_entry.get().strip()

    if not compound_name:
        messagebox.showerror("Error", "No compound name selected, please enter a compound name")
        return

    try:
        matches = library.find_by_name(compound_name)
        print("Matches:", matches)

        if not matches:
            result_text.delete("1.0", tk.END)

            result_text.insert(tk.END, f"No Matches found for:\n{compound_name}\n")
            return

        #First matching pairs

        feature_id = matches[0]

        print("Selected feature ID:", feature_id)

        feature = library.get_feature(feature_id)

        print("Feature:", feature)

        #Get Properties
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, "Compound Information\n")
        result_text.insert(tk.END, f"Compound: {compound_name}\n")
        result_text.insert(tk.END, f"Feature: {feature}\n")

        #Display Library properties
        result_text.insert(tk.END, f"Name: {library.names.get(feature_id)}\n")
        result_text.insert(tk.END, f"Molecular Mass: {library.molecular_mass.get(feature_id)}\n")
        result_text.insert(tk.END, f"m/z: {library.mzs.get(feature_id)}\n")
        result_text.insert(tk.END, f"Metaboscape formulas: {library.formula_metaboscape.get(feature_id)}")
        result_text.insert(tk.END, f"SIRIUS formulas: {library.formula_sirius.get(feature_id)}")
        result_text.insert(tk.END, f"SMILES: {library.smiles.get(feature_id)}\n")
        result_text.insert(tk.END, f"InChI: {library.inchis.get(feature_id)}\n")
        result_text.insert(tk.END, f"Retention time: {library.retention_times_in_seconds.get(feature_id)}\n")
        result_text.insert(tk.END, f"Collisional cross section: {library.collisional_cross_sections.get(feature_id)}")

    except Exception as e:
        messagebox.showerror("Error", str(e))

#GUI
root = tk.Tk()
root.geometry("800x700")
root.title("Library Info")
root.columnconfigure(0, weight=1)

#Database section
database_frame = tk.LabelFrame(root, text="Database", padx=10, pady=10)
database_frame.grid(row=0, column=0, padx=20, pady=20, sticky="ew")
database_frame.columnconfigure(1, weight=1)

#SQLite file
tk.Label(database_frame, text="SQLite path:").grid(row=0, column=0, padx=10, pady=5, sticky="w")
db_entry = tk.Entry(database_frame)
db_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")

#Brwose button
tk.Button(database_frame, text="Browse", command=select_database).grid(row=0, column=2, padx=10, pady=5, sticky="w")

#Load Button
tk.Button(database_frame, text="Load Database", command=load_database).grid(row=1, column=1, padx=10, pady=10, sticky="w")

#Status
status_label = tk.Label(database_frame, text="No database loaded.")
status_label.grid(row=1, column=2, padx=10, pady=10)

#Compound Section
compound_frame = tk.LabelFrame(root, text="Compound", padx=10, pady=10)
compound_frame.grid(row=1, column=0, padx=20, pady=10, sticky="ew")
tk.Label(compound_frame, text="Compound name:").grid(row=0, column=0, padx=10, pady=5, sticky="w")

#Compound name entry
compound_entry = tk.Entry(compound_frame)
compound_entry.grid(row=0, column=1, padx=10, pady=5, sticky="ew")

#Search button
tk.Button(compound_frame, text="Search", command=search_compound).grid(row=0, column=2, padx=10, pady=5)

#Properties button
tk.Button(compound_frame, text="Show Properties", command=show_compound_properties).grid(row=1, column=1, padx=10, pady=10, sticky="w")

#Result Section
result_frame = tk.LabelFrame(root, text="Compound Properties", padx=10, pady=10)
result_frame.grid(row=2, column=0, padx=20, pady=20, sticky="nsew")
root.rowconfigure(2, weight=1)
result_frame.rowconfigure(0, weight=1)

#Box to display info
result_text=tk.Text(result_frame, wrap="word")
result_text.grid(row=0, column=0, sticky="nsew")

#Scroll bar
scroll_bar = tk.Scrollbar(result_frame, orient="vertical", command=result_text.yview)
scroll_bar.grid(row=0, column=1, sticky="nse")
result_text.config(yscrollcommand=scroll_bar.set)

#Library info button
tk.Button(root, text="Show Library information", command=display_library_info).grid(row=3, column=0, pady=15)

#Start GUI
root.mainloop()