import tkinter as tk
from tkinter import filedialog, messagebox
from msIO.feature_managers.db import Library
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  #Got to know about it from Matplotlib documentation
from matplotlib import pyplot as plt


class GuiManager:
    library = None
    canvas = None
    db_file = None

    #Selecting the SQLite file
    def select_database(self):
        path = filedialog.askopenfilename(title="Select SQLite database",filetypes=[("sqlite files", "*.sqlite"), ("SQLite databse", "*.db"), ("All files", "*.*")])

        if not path:
            return
        try:
            self.library = Library(path)
            database_label.config(text=f"Database: {path}")
            feature_label.config(text=f"Number of features: {len(self.library.feature_ids)}")
            messagebox.showinfo("Success", "SQLite library opened successfully!")
        except Exception as e:
            messagebox.showerror("Error", f"Could not open databse:\n\n{e}")

    def plot_compound(self):
        if self.library is None:
            messagebox.showwarning("No database", "Please select a database first in SQLite format")
            return
        if len(self.library.feature_ids) == 0:
            messagebox.showwarning("Empty library", "There are no features")
            return
        try:
            feature_id = self.library.feature_ids[0]
            print("selected features ID: ", feature_id)

            fig, axs = plt.subplots(1, 2, figsize=(10, 5))

            if self.canvas is not None:
                self.canvas.get_tk_widget().destroy()

            self.canvas = FigureCanvasTkAgg(fig, master=plot_frame)

            self.canvas.draw()
            self.canvas.get_tk_widget().pack(fill = tk.BOTH, expand = True)
            print("Plot displayed successfully")

        except Exception as e:
            messagebox.showerror("Error", f"Could not plot:\n\n{e}")
            print("Plot error:", e)

mgr = GuiManager()

#Main window
root = tk.Tk()
root.title("msIO Library viewer")
root.geometry("1100x700")

#Top frame
top_frame = tk.Frame(root)
top_frame.pack(fill = tk.X, padx = 10, pady = 10)

#Database button
select_button = tk.Button(top_frame, text="Select database", command=mgr.select_database)
select_button.pack(side = tk.LEFT, padx=5)

#PLOT Button
plot_button = tk.Button(top_frame, text="Plot", command=mgr.plot_compound)
plot_button.pack(side = tk.LEFT, padx=5)

#Database Info
database_label = tk.Label(root, text="Database: None", anchor="w")
database_label.pack(fill = tk.X, padx = 15)
feature_label = tk.Label(root, text="Number of features: -", anchor="w")
feature_label.pack(fill = tk.X, padx = 15, pady = (0,10))

#Frame for the plot
plot_frame = tk.Frame(root, bd=2, relief=tk.SUNKEN)
plot_frame.pack(fill = tk.BOTH, expand=True, padx = 10, pady = 10)

#Start
root.mainloop()