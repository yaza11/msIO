# msIO

msIO reads mass spectrometry data from several tools, combines it into one SQLite database, and provides a GUI for matching measured features against spectral libraries and annotating them.

## Features

- **Readers** for MS data formats and tool exports:
  - MGF (MS1/MS2 spectra)
  - MSP spectral libraries
  - MetaboScape feature tables and `.mca` annotation files
  - MZmine feature tables and `.graphml` networks
  - GNPS molecular networks
  - SIRIUS summaries (formula, compound and CANOPUS identifications)
- **SQLite storage**: combine a feature table with MGF, GNPS and SIRIUS results into one database. MSP libraries can be turned into library databases too.
- **Annotation GUI**: match features against one or more libraries (precursor m/z, MS2 similarity), filter the results, look at MS1/MS2 spectra and structures, and annotate features, peaks, adducts, isomers and series.

## Installation

Python 3.10 or newer is needed. msIO needs [LipidCalculator](https://github.com/yaza11/LipidCalculator), which is not on PyPI yet and has to be installed first:

```bash
pip install git+https://github.com/yaza11/LipidCalculator.git@v0.3.7
pip install msIO
```

The latest development version can be installed from GitHub with `pip install git+https://github.com/yaza11/msIO.git` instead.

For Windows, there is also a standalone GUI that needs no Python installation: download `msIO_GUI-windows.zip` from the [releases](https://github.com/yaza11/msIO/releases), extract it and run `msIO_GUI.exe`.

## Graphical user interface

Start the GUI with:

```bash
python -m msIO.gui.main
```

The window has these tabs:

| Tab | Purpose |
| --- | --- |
| MZmine/MetaboScape → SQLite | Combine a feature table with MGF, GNPS and SIRIUS output into a measured-feature database |
| MSP library → SQLite | Convert one or more `.msp` libraries into a library database |
| Library viewer | Search a library database by compound name and plot entries |
| Annotation | Match measured features against libraries and annotate them |
| Guide | A suggested strategy for making annotations |

You can also run each tab on its own by running its module in `msIO/gui/`.

[docs/GUI_features.md](docs/GUI_features.md) explains the annotation GUI in detail, with screenshots.

### Standalone executable

Ready-built Windows executables are attached to the [releases](https://github.com/yaza11/msIO/releases). To build one yourself (no Python needed to run it), install PyInstaller in the environment where msIO and its dependencies are installed and run:

```bash
pyinstaller --noconfirm --clean LipidAnnotations.spec
```

The build goes to `dist/msIO_GUI/`. The GitHub workflow builds the releases from the same spec file.

## Python usage

### Build a measured-feature database

```python
from msIO import MgfImportManager
from msIO.feature_managers.combined import ProjectImportManager
from msIO.feature_managers.gnps import GnpsImportManager
from msIO.feature_managers.metaboscape import MetaboscapeImportManager
from msIO.feature_managers.sirius import SiriusImportManager
from msIO.sql.session import initiate_db

manager = ProjectImportManager(
    metaboscape_manager=MetaboscapeImportManager("height.csv"),
    mgf_manager=MgfImportManager("ms1_ms2.mgf", accumulate_spectra=False),
    gnps_manager=GnpsImportManager("gnps"),
    sirius_manager=SiriusImportManager(path_folder_export="sirius"),
)

initiate_db("features.sqlite")
manager.to_sql("features.sqlite", feature_ids=manager.feature_ids)
```

### Build a library database from MSP files

```python
from msIO.sql.create_library import write_lib_from_msp_files

write_lib_from_msp_files("library.sqlite", ["lib1.msp", "lib2.msp"])
```

### Query databases

```python
from msIO.feature_managers.db import FeatureManagerDB, Library

features = FeatureManagerDB("features.sqlite")
spectrum = features.get_ms_spectrum(features.feature_ids[0], level=2)

lib = Library("library.sqlite")
entry_id = lib.find_by_name("OH-GTGT-0a")[0]
lib.plot_compound_overview(entry_id)
```

## Package layout

| Module | Contents |
| --- | --- |
| `msIO.list_of_ions` | Readers for peak lists: MGF, MSP, MCA and library spectra |
| `msIO.features` | Feature classes for each data source |
| `msIO.feature_managers` | Import managers for each source, plus database access (`FeatureManagerDB`, `Library`) |
| `msIO.sql` | SQLAlchemy schema, sessions and library creation |
| `msIO.annotations` | Annotation handling |
| `msIO.gui` | Tkinter GUIs |

## License

See [LICENSE](LICENSE).
