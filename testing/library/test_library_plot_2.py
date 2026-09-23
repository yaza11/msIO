import tkinter as tk
from tkinter import filedialog, messagebox
from msIO.feature_managers.db import Library
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

class LibraryGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("msIO Library Viewer")
        self.root.geometry("1200x800")

        self.library = None
        self.search_results = []
        self.canvas = None

        #Creating GUI

        self.create_gui()

    def create_gui(self):
        database_frame = tk.Frame(self.root)
        database_frame.pack(fill = tk.X, padx=10, pady = 10)
        self.select_database_button = tk.Button(database_frame, text = "Select SQLite Database", command = self.select_database)
        self.select_database_button.pack(side = tk.LEFT, padx=5)
        self.database_label = tk.Label(database_frame, text = "Database:None", anchor = "w")
        self.database_label.pack(side = tk.LEFT, padx=10)

        #Feature info
        self.feature_label = tk.Label(self.root, text="Number of features: -", anchor = "w")
        self.feature_label.pack(side = tk.LEFT, padx=15)

        #Compound search
        search_frame = tk.LabelFrame(self.root, text = "Compound Search")
        search_frame.pack(fill = tk.X, padx=10, pady = 10)

        name_label = tk.Label(search_frame, text = "Name:")
        name_label.grid(row = 0, column = 0, padx=5, pady=10)

        self.name_entry = tk.Entry(search_frame, width=40)
        self.name_entry.grid(row = 0, column = 1, padx=5, pady=10)

        #Search Button
        self.search_button = tk.Button(search_frame, text = "Search", command = self.search_compound)
        self.search_button.grid(row = 0, column = 2, padx=5, pady=10)

        #Search Results
        result_label = tk.Label(search_frame, text = "Result:")
        result_label.grid(row = 1, column = 0, padx=5, pady=5, sticky = "nw")
        self.result_listbox = tk.Listbox(search_frame, height = 5, width = 80)
        self.result_listbox.grid(row = 1, column = 1, columnspan = 2, padx=5, pady=5, sticky = "ew")

        #Plot Button
        self.plot_button = tk.Button(search_frame, text = "Plot selected compound", command = self.plot_selected_compound)
        self.plot_button.grid(row = 2, column = 1, padx=5, pady=10, sticky = "w")

        #Plot Area
        self.plot_frame = tk.LabelFrame(self.root, text="Compound Overview")
        self.plot_frame.pack(fill = tk.BOTH, expand=True, padx=10, pady=10)

    #Select database
    def select_database(self):
        db_file = filedialog.askopenfilename(title="Select SQLite library", filetypes=[("SQLite files", "*.sqlite"), ("SQLite database", "*.db"), ("All files", "*")])

        if not db_file:
            return
        try:
            self.library = Library(db_file)
            #Display database path
            self.database_label.config(text = f"Database: {db_file}")

            #Display features
            self.feature_label.config(text=f"Number of features:" f"{len(self.library.feature_ids)}")

            #Remove previous search results
            self.result_listbox.delete(0, tk.END)

            #Reset stored search result
            self.search_results = []
            messagebox.showinfo("Success", "Successfully opened SQLite library")
            print("Library opened successfully")
            print("Number of features: ", len(self.library.feature_ids))

        except Exception as e:
            messagebox.showerror("Error", f"Could not open the database:\n\n{e}")
            print("Database error:", e)


    def search_compound(self):
        if self.library is None:
            messagebox.showerror("No database", "Please select a SQLite database")
            return

        #Get compound name
        compound_name = self.name_entry.get().strip()
        if not compound_name:
            messagebox.showerror("Error", "Please enter a compound name")
            return
        try:
            print()
            print("Searching for:", compound_name)

            #Using find_by_name as suggested by Yannick
            self.search_results: list[int] = self.library.find_by_name(compound_name)
            print("Search Result:")
            print(self.search_results)

            #Clear old result
            self.result_listbox.delete(0, tk.END)

            #Checking if nothing is found
            if not self.search_results:
                messagebox.showerror("No results", f"No compound found for:\n{compound_name}")
                return

            #Display results
            for result in self.search_results:
                self.result_listbox.insert(tk.END, str(result))
            print("Number of results:", len(self.search_results))

        except Exception as e:
            messagebox.showerror("Error", f"Could not search for compound:\n\n{e}")
            print("Search error:", e)

    def plot_selected_compound(self):
        if self.library is None:
            messagebox.showerror("No database", "Please select a SQLite database")
            return

        if not self.search_results:
            messagebox.showerror("No results", "Please search for a compound first")
            return

        #Find selected item
        selected_index: int = self.result_listbox.curselection()[0]
        print(selected_index)
        if selected_index is None:
            messagebox.showerror("No results", "Please select a compound first")
            return

        selected_result = self.search_results[selected_index]
        print()
        print("Selected result:")
        print(selected_result)

        try:
            #Using selected result as feature ID
            feature_id = selected_result
            print("Feature ID used for plotting:", feature_id)

            # Remove the previous canvas
            if self.canvas is not None:
                self.canvas.get_tk_widget().destroy()
                self.canvas = None

            #Creating matplotlib figure
            fig, axs = plt.subplots(1, 2, figsize = (10, 5))
            plt.close(fig)

            #Asking msIO to create compound overview
            self.library.plot_compound_overview(feature_id, axs=axs)

            #Putting Matplotlib inside the Tkinter
            self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
            self.canvas.draw()
            self.canvas.get_tk_widget().pack(fill = tk.BOTH, expand=True)
            print("Compound overview displayed")

        except Exception as e:
            plt.close("all")
            messagebox.showerror("Plot Error", "Could not plot selected compound")

            print("Plot error:", e)
            print(e)

#Start Program

root = tk.Tk()
app = LibraryGUI(root)
root.mainloop()
