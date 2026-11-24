"""
msIO main window: all GUIs of msIO.gui as tabs of one window.

    Annotation              match measured features against spectral libraries and annotate them (msIO.gui.annotation)
    Library viewer          search a library database by compound name and plot entries (library_search_by_name)
    Feature table → SQLite  combine a MetaboScape / MZmine feature table with MGF, GNPS and SIRIUS (csv_to_sqlite_gui)
    MSP → SQLite            convert msp spectral libraries into library databases (msp_to_sqlite_gui)
    Guide                   suggested strategy for making annotations (text file "about")

Each tab can also be started on its own by running its module.
"""
import os
import sys

# a windowed executable (PyInstaller, console=False) has no stdout / stderr, printing would fail
if sys.stdout is None:
    sys.stdout = open(os.path.join(os.path.expanduser("~"), "msio_converter.log"), "a", buffering=1)
if sys.stderr is None:
    sys.stderr = open(os.path.join(os.path.expanduser("~"), "msio_converter_errors.log"), "a", buffering=1)

import tkinter as tk
from tkinter import ttk

from msIO import __version__
from msIO.gui.annotation import AnnotationGUI
from msIO.gui.csv_to_sqlite_gui import ConverterApp
from msIO.gui.library_search_by_name import LibraryNameSearch
from msIO.gui.msp_to_sqlite_gui import MspConverterApp

PATH_GUIDE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "about")


class MainWindow:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title(f"msIO {__version__}")
        # 80 % of the screen, the annotation plots need the space
        width: int = int(root.winfo_screenwidth() * 0.8)
        height: int = int(root.winfo_screenheight() * 0.8)
        root.geometry(f"{width}x{height}+10+10")
        root.protocol("WM_DELETE_WINDOW", self.on_close)

        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill="both", expand=True)

        self.feature_table_converter = ConverterApp(self._add_form_tab("MZmine/MetaboScape --> SQLite"))
        self.msp_converter = MspConverterApp(self._add_form_tab("MSP library --> SQLite"))
        self.library_viewer = LibraryNameSearch(self._add_tab("Library viewer"))
        self.annotation = AnnotationGUI(self._add_tab("Annotation"))
        if os.path.isfile(PATH_GUIDE):
            self._add_guide(self._add_tab("Guide"))

    def _add_tab(self, title: str) -> ttk.Frame:
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text=title)
        return frame

    def _add_form_tab(self, title: str, width: int = 900) -> ttk.Frame:
        """Tab with a centered column of fixed width, forms stretched over the whole screen are hard to read."""
        column = ttk.Frame(self._add_tab(title))
        column.pack(anchor="n", pady=10)
        column.columnconfigure(0, minsize=width)  # the forms are gridded in column 0, the height follows the content
        return column

    @staticmethod
    def _add_guide(frame: ttk.Frame) -> None:
        with open(PATH_GUIDE, encoding="utf-8") as f:
            content: str = f.read()
        text = tk.Text(frame, wrap="word", padx=15, pady=15, relief="flat")
        scrollbar = ttk.Scrollbar(frame, command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)
        text.insert("1.0", content)
        text.configure(state="disabled")

    def on_close(self) -> None:
        """Ask before discarding annotation changes that were not written to the database."""
        if self.annotation.can_close():
            self.root.destroy()


def main() -> None:
    root = tk.Tk()
    MainWindow(root)
    root.mainloop()


if __name__ == "__main__":
    main()
