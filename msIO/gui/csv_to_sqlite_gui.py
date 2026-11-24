"""GUI for combining a feature table (MetaboScape or MZmine) with MGF, GNPS and SIRIUS
results into a SQLite database.

Originally written by Avin Jain.
"""
# noinspection PyUnresolvedReferences
import msIO.gui.executable  # noqa: F401  (must stay first)
import glob
import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox

from msIO import MgfImportManager
from msIO.feature_managers.combined import ProjectImportManager
from msIO.feature_managers.gnps import GnpsImportManager
from msIO.feature_managers.metaboscape import MetaboscapeImportManager
from msIO.feature_managers.sirius import SiriusImportManager
from msIO.sql.session import initiate_db

SOURCE_METABOSCAPE = 'MetaboScape'
SOURCE_MZMINE = 'MZmine'

# file name patterns of a default MZmine export, the first match is used
MZMINE_EXPORT_PATTERNS = {
    'feature_table': ['*_full_feature_table.csv'],
    'mgf': ['*_sirius.mgf', '*.mgf'],
    'gnps': ['*_networking_iimn.graphml', '*_networking_fbmn.graphml', '*.graphml'],
}

PAD = dict(padx=10, pady=5)


class PathRow:
    """Label, entry and browse button in one row of a frame."""

    def __init__(self, frame: tk.Widget, row: int, label: str, browse):
        tk.Label(frame, text=label).grid(row=row, column=0, sticky='w', **PAD)
        self.entry = tk.Entry(frame)
        self.entry.grid(row=row, column=1, sticky='ew', **PAD)
        tk.Button(frame, text='Browse...', command=lambda: self.set(browse())).grid(row=row, column=2, **PAD)

    def get(self) -> str:
        return self.entry.get().strip()

    def set(self, path: str) -> None:
        if path:
            self.entry.delete(0, tk.END)
            self.entry.insert(0, os.path.normpath(path))


def _section(root: tk.Tk, row: int, title: str) -> tk.LabelFrame:
    frame = tk.LabelFrame(root, text=title, padx=10, pady=10)
    frame.grid(row=row, column=0, padx=20, pady=10, sticky='ew')
    frame.columnconfigure(1, weight=1)
    return frame


def _find_first(folder: str, patterns: list[str]) -> str | None:
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.join(glob.escape(folder), pattern)))
        if matches:
            return matches[0]
    return None


class ConverterApp:
    def __init__(self, root: tk.Misc):
        """root: a window or a frame to embed the converter in."""
        self.root = root
        if isinstance(root, (tk.Tk, tk.Toplevel)):
            root.title('Feature table to SQLite converter')
        root.columnconfigure(0, weight=1)

        # feature table
        frame = _section(root, 0, 'Feature table')
        self.source = tk.StringVar(value=SOURCE_MZMINE)
        source_frame = tk.Frame(frame)
        source_frame.grid(row=0, column=1, sticky='w')
        for source in (SOURCE_MZMINE, SOURCE_METABOSCAPE):
            tk.Radiobutton(source_frame, text=source, value=source, variable=self.source).pack(side='left', padx=(0, 10))
        tk.Button(frame, text='Fill from MZmine export folder...', command=self.fill_from_mzmine_folder) \
            .grid(row=1, column=1, sticky='w', **PAD)
        self.feature_table = PathRow(frame, 2, 'CSV file', lambda: self._ask_file('CSV', '*.csv'))
        self.intensity = tk.StringVar(value='area')
        intensity_frame = tk.Frame(frame)
        intensity_frame.grid(row=3, column=1, sticky='w')
        tk.Label(intensity_frame, text='MZmine intensity:').pack(side='left')
        for value in ('area', 'height'):
            tk.Radiobutton(intensity_frame, text=value, value=value, variable=self.intensity).pack(side='left')

        # MGF
        frame = _section(root, 1, 'MGF (optional)')
        self.mgf = PathRow(frame, 0, 'MGF file', lambda: self._ask_file('MGF', '*.mgf'))
        self.accumulate = tk.BooleanVar(value=True)
        tk.Checkbutton(frame, text='Accumulate spectra', variable=self.accumulate,
                       command=self._toggle_mz_resolution).grid(row=1, column=1, sticky='w', **PAD)
        tk.Label(frame, text='m/z resolution').grid(row=2, column=0, sticky='w', **PAD)
        self.mz_resolution = tk.Entry(frame)
        self.mz_resolution.insert(0, '40000')
        self.mz_resolution.grid(row=2, column=1, sticky='ew', **PAD)

        # GNPS
        frame = _section(root, 2, 'Molecular network (optional)')
        self.gnps = PathRow(frame, 0, 'GNPS folder or\nMZmine .graphml', self._ask_gnps)

        # SIRIUS
        frame = _section(root, 3, 'SIRIUS (optional)')
        self.sirius = PathRow(frame, 0, 'Export folder', lambda: filedialog.askdirectory(title='Select SIRIUS folder'))
        tk.Label(frame, text='Export tag').grid(row=1, column=0, sticky='w', **PAD)
        self.sirius_tag = tk.Entry(frame)
        self.sirius_tag.grid(row=1, column=1, sticky='ew', **PAD)

        # output
        frame = _section(root, 4, 'Output')
        self.sqlite = PathRow(frame, 0, 'SQLite file', lambda: filedialog.asksaveasfilename(
            title='Select SQLite database', defaultextension='.sqlite',
            filetypes=[('SQLite files', '*.sqlite *.db'), ('All files', '*.*')]))

        self.convert_button = tk.Button(root, text='Convert', width=15, command=self.convert)
        self.convert_button.grid(row=5, column=0, pady=(10, 0))
        self.status = tk.Label(root, text='', fg='gray30')
        self.status.grid(row=6, column=0, pady=(5, 15))

    # --- browsing -------------------------------------------------------------------------
    @staticmethod
    def _ask_file(name: str, pattern: str) -> str:
        return filedialog.askopenfilename(
            title=f'Select {name} file', filetypes=[(f'{name} files', pattern), ('All files', '*.*')])

    def _ask_gnps(self) -> str:
        if self.source.get() == SOURCE_MZMINE:
            return self._ask_file('GraphML', '*.graphml')
        return filedialog.askdirectory(title='Select GNPS folder')

    def _toggle_mz_resolution(self) -> None:
        self.mz_resolution.config(state='normal' if self.accumulate.get() else 'disabled')

    def fill_from_mzmine_folder(self) -> None:
        folder = filedialog.askdirectory(title='Select MZmine export folder')
        if not folder:
            return
        self.source.set(SOURCE_MZMINE)
        found = {key: _find_first(folder, patterns) for key, patterns in MZMINE_EXPORT_PATTERNS.items()}
        if found['feature_table'] is None:
            messagebox.showwarning('Not found', f'No *_full_feature_table.csv found in\n{folder}')
        self.feature_table.set(found['feature_table'])
        self.mgf.set(found['mgf'])
        self.gnps.set(found['gnps'])
        if not self.sqlite.get() and found['feature_table']:
            self.sqlite.set(found['feature_table'].removesuffix('_full_feature_table.csv') + '.sqlite')

    # --- conversion -----------------------------------------------------------------------
    def _collect_settings(self) -> dict:
        """Validate the inputs, raise ValueError with a user-readable message if something is wrong."""
        settings = dict(
            source=self.source.get(),
            feature_table=self.feature_table.get(),
            intensity=self.intensity.get(),
            mgf=self.mgf.get(),
            mz_resolution=None,
            gnps=self.gnps.get(),
            sirius=self.sirius.get(),
            sirius_tag=self.sirius_tag.get() or None,
            sqlite=self.sqlite.get(),
        )
        if not settings['feature_table']:
            raise ValueError('Please select a feature table (CSV file).')
        if not settings['sqlite']:
            raise ValueError('Please select an output SQLite file.')
        for key in ('feature_table', 'mgf', 'gnps', 'sirius'):
            if settings[key] and not os.path.exists(settings[key]):
                raise ValueError(f'Path does not exist:\n{settings[key]}')
        if settings['mgf'] and self.accumulate.get():
            try:
                settings['mz_resolution'] = int(self.mz_resolution.get())
            except ValueError:
                raise ValueError('The m/z resolution must be an integer, e.g. 40000.')
        return settings

    def convert(self) -> None:
        try:
            settings = self._collect_settings()
        except ValueError as e:
            messagebox.showerror('Error', str(e))
            return
        if os.path.exists(settings['sqlite']) and not messagebox.askyesno(
                'Overwrite?', f'{settings["sqlite"]}\nalready exists. Overwrite it?'):
            return

        self.convert_button.config(state='disabled')
        threading.Thread(target=self._run, args=(settings,), daemon=True).start()

    def _set_status(self, text: str) -> None:
        self.root.after(0, lambda: self.status.config(text=text))

    def _run(self, settings: dict) -> None:
        try:
            run_conversion(settings, self._set_status)
        except Exception as e:  # report any failure in the GUI
            error = f'{type(e).__name__}: {e}'
            self.root.after(0, lambda: self._finish('Conversion failed.', error))
        else:
            self.root.after(0, lambda: self._finish('Done.'))

    def _finish(self, status: str, error: str = None) -> None:
        self.status.config(text=status)
        self.convert_button.config(state='normal')
        if error is None:
            messagebox.showinfo('Success', 'Database written successfully.')
        else:
            messagebox.showerror('Error', error)


def run_conversion(settings: dict, report=print) -> None:
    """Read all selected inputs and write them to the SQLite file in settings['sqlite']."""
    report('Reading feature table...')
    if settings['source'] == SOURCE_MZMINE:
        feature_table = MetaboscapeImportManager.from_mzmine(settings['feature_table'], settings['intensity'])
    else:
        feature_table = MetaboscapeImportManager(settings['feature_table'])

    mgf = None
    if settings['mgf']:
        report('Reading MGF file...')
        mgf = MgfImportManager(
            settings['mgf'],
            accumulate_spectra=settings['mz_resolution'] is not None,
            mz_resolution=settings['mz_resolution'],
        )

    gnps = None
    if path_gnps := settings['gnps']:
        report('Reading molecular network...')
        if os.path.isdir(path_gnps):
            gnps = GnpsImportManager(path_gnps_folder=path_gnps)
        elif settings['source'] == SOURCE_MZMINE:
            gnps = GnpsImportManager.from_mzmine(path_gnps)
        else:
            gnps = GnpsImportManager(path_file_gnps_graphml=path_gnps)

    sirius = None
    if settings['sirius']:
        report('Reading SIRIUS results...')
        sirius = SiriusImportManager(path_folder_export=settings['sirius'], export_tag=settings['sirius_tag'])

    project = ProjectImportManager(
        metaboscape_manager=feature_table,
        mgf_manager=mgf,
        gnps_manager=gnps,
        sirius_manager=sirius,
    )

    report(f'Writing {len(feature_table.feature_ids)} features to database...')
    initiate_db(settings['sqlite'])
    project.to_sql(settings['sqlite'], feature_ids=feature_table.feature_ids)


def main() -> None:
    root = tk.Tk()
    ConverterApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
