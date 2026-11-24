"""GUI for converting one or more msp spectral libraries into a SQLite library database."""
# noinspection PyUnresolvedReferences
import msIO.gui.executable  # noqa: F401  (must stay first)
import os
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from msIO.sql.create_library import write_lib_from_msp_files

PAD = dict(padx=10, pady=5)


class MspConverterApp:
    def __init__(self, root: tk.Misc):
        """root: a window or a frame to embed the converter in."""
        self.root = root
        if isinstance(root, (tk.Tk, tk.Toplevel)):
            root.title('MSP library to SQLite converter')
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

        # input files
        frame = tk.LabelFrame(root, text='MSP libraries', padx=10, pady=10)
        frame.grid(row=0, column=0, padx=20, pady=(20, 10), sticky='nsew')
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        self.files = tk.Listbox(frame, height=6, selectmode=tk.EXTENDED)
        self.files.grid(row=0, column=0, rowspan=2, sticky='nsew', **PAD)
        tk.Button(frame, text='Add...', width=10, command=self.add_files).grid(row=0, column=1, sticky='n', **PAD)
        tk.Button(frame, text='Remove', width=10, command=self.remove_files).grid(row=1, column=1, sticky='n', **PAD)

        # output
        frame = tk.LabelFrame(root, text='Output', padx=10, pady=10)
        frame.grid(row=1, column=0, padx=20, pady=10, sticky='ew')
        frame.columnconfigure(1, weight=1)
        tk.Label(frame, text='SQLite file').grid(row=0, column=0, sticky='w', **PAD)
        self.sqlite = tk.Entry(frame)
        self.sqlite.grid(row=0, column=1, sticky='ew', **PAD)
        tk.Button(frame, text='Browse...', command=self.select_output).grid(row=0, column=2, **PAD)

        self.convert_button = tk.Button(root, text='Convert', width=15, command=self.convert)
        self.convert_button.grid(row=2, column=0, pady=(10, 0))
        self.progress = ttk.Progressbar(root, maximum=1.)
        self.progress.grid(row=3, column=0, padx=20, pady=(10, 0), sticky='ew')
        self.status = tk.Label(root, text='', fg='gray30')
        self.status.grid(row=4, column=0, pady=(5, 15))

    def _set_output(self, path: str) -> None:
        self.sqlite.delete(0, tk.END)
        self.sqlite.insert(0, os.path.normpath(path))

    def add_files(self) -> None:
        paths = filedialog.askopenfilenames(
            title='Select MSP files', filetypes=[('MSP files', '*.msp'), ('All files', '*.*')])
        current = set(self.files.get(0, tk.END))
        for path in map(os.path.normpath, paths):
            if path not in current:
                self.files.insert(tk.END, path)
        if paths and not self.sqlite.get():
            self._set_output(os.path.splitext(paths[0])[0] + '.sqlite')

    def remove_files(self) -> None:
        for idx in reversed(self.files.curselection()):
            self.files.delete(idx)

    def select_output(self) -> None:
        path = filedialog.asksaveasfilename(
            title='Select SQLite database', defaultextension='.sqlite',
            filetypes=[('SQLite files', '*.sqlite *.db'), ('All files', '*.*')])
        if path:
            self._set_output(path)

    def convert(self) -> None:
        msp_files = list(self.files.get(0, tk.END))
        db_file = self.sqlite.get().strip()
        if not msp_files:
            messagebox.showerror('Error', 'Please add at least one MSP file.')
            return
        if not db_file:
            messagebox.showerror('Error', 'Please select an output SQLite file.')
            return
        if missing := [f for f in msp_files if not os.path.isfile(f)]:
            messagebox.showerror('Error', 'File(s) not found:\n' + '\n'.join(missing))
            return
        if os.path.abspath(db_file) in map(os.path.abspath, msp_files):
            messagebox.showerror('Error', 'The output file must not be one of the input files.')
            return
        if os.path.exists(db_file) and not messagebox.askyesno(
                'Overwrite?', f'{db_file}\nalready exists. Overwrite it?'):
            return

        self.convert_button.config(state='disabled')
        threading.Thread(target=self._run, args=(db_file, msp_files), daemon=True).start()

    def _run(self, db_file: str, msp_files: list[str]) -> None:
        t_start = time.time()
        n_files = len(msp_files)
        # t_last: time of the last progress report, a new file starts when the previous one is reported as done
        state = {'file': 0, 'name': None, 't_file': t_start, 't_last': t_start}

        def progress(lib_name: str, n_entries: int, fraction: float) -> None:
            if lib_name != state['name']:
                state.update(file=state['file'] + 1, name=lib_name, t_file=state['t_last'])
            state['t_last'] = time.time()
            text = f'[{state["file"]}/{n_files}] {lib_name}: {n_entries:,} entries'
            elapsed = time.time() - state['t_file']
            if 0 < fraction < 1:
                text += f' (~{elapsed / fraction * (1 - fraction):.0f} s left)'
            self.root.after(0, self._update_progress, text, fraction)

        try:
            n_entries = write_lib_from_msp_files(db_file, msp_files, progress=progress)
        except Exception as e:  # report any failure in the GUI
            self.root.after(0, self._finish, f'Conversion failed: {type(e).__name__}: {e}')
        else:
            self.root.after(0, self._finish, None, time.time() - t_start, n_entries)

    def _update_progress(self, text: str, fraction: float) -> None:
        self.status.config(text=text)
        self.progress.config(value=fraction)

    def _finish(self, error: str | None, seconds: float = 0., n_entries: int = 0) -> None:
        self.convert_button.config(state='normal')
        if error is not None:
            self.status.config(text='Conversion failed.')
            messagebox.showerror('Error', error)
        else:
            self.status.config(text=f'Wrote {n_entries:,} entries in {seconds:.0f} s.')
            messagebox.showinfo('Success', 'Library database written successfully.')


def main() -> None:
    root = tk.Tk()
    root.geometry('700x400')
    MspConverterApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
