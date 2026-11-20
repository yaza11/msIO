"""Write msp spectral libraries into a SQLite library database.

Entries are parsed from the msp file as a stream and inserted in bulk with SQLAlchemy core, which is much faster
than building ORM objects for every entry and peak. The rows are the same as those of the ORM objects created by
MSPReader.create_feature (one FeatureCombined per entry with a MetaboScape, MGF and SIRIUS part).
"""
import os
import re
from typing import Callable, Iterator

from sqlalchemy import insert

from msIO.features.combined import FeatureCombined
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.mgf import FeatureMgf, MsSpec
from msIO.features.sirius import FeatureSirius, CompoundCandidate
from msIO.list_of_ions.base import PeakList, PeakFeature
from msIO.list_of_ions.read_msp import _parse_lines
from msIO.sql.session import initiate_db, get_engine

# models in the order they have to be inserted (referenced tables first)
_MODELS = (FeatureCombined, PeakList, PeakFeature, FeatureMgf, MsSpec, FeatureMetaboScape, FeatureSirius,
           CompoundCandidate)

_PEAK_SEPARATOR = re.compile(r'[\t ,]+')


def _value_or_none(val):
    """Same as in MSPReader.create_feature: negative numbers and NaN are treated as missing."""
    if isinstance(val, int | float) and not (val >= 0):
        return None
    return val


def _iter_entries(path_file: str) -> Iterator[tuple[list[str], list[str], int]]:
    """Yield (metadata lines, peak lines, bytes read so far) for every entry in the msp file."""
    meta, peaks = [], []
    n_bytes = 0
    with open(path_file, 'rb') as f:
        for raw in f:
            n_bytes += len(raw)
            line = raw.decode('utf-8', errors='replace').strip('\r\n')
            if not line.strip():  # entries are separated by blank lines
                if meta or peaks:
                    yield meta, peaks, n_bytes
                    meta, peaks = [], []
            elif line[0].isdigit():
                peaks.append(line)
            else:
                meta.append(line)
    if meta or peaks:
        yield meta, peaks, n_bytes


def _parse_peak(line: str) -> tuple[float, float, str | None]:
    """'mz intensity "annotation"' (tab, comma or space separated, annotation optional)"""
    annotation = None
    if '"' in line:
        line, annotation = line.split('"', 1)
        annotation = annotation.rstrip('"')
    mz, intensity = _PEAK_SEPARATOR.split(line.strip())[:2]
    return float(mz), float(intensity), annotation


def write_lib_from_msp_files(
        db_file: str,
        msp_files: list[str],
        commit_at_latest_after: int = 10_000,
        progress: Callable[[str, int, float], None] = None
) -> int:
    """Write all entries of the msp files into a new library database at db_file (an existing file is replaced).

    commit_at_latest_after: number of entries that are inserted at once
    progress: if given, called as progress(library_name, n_entries_written, fraction_of_file_read) for every batch

    Returns the number of entries written.
    """
    initiate_db(db_file)
    engine = get_engine(db_file)

    rows: dict[type, list[dict]] = {model: [] for model in _MODELS}
    f_id = 0  # running index for features, also used as id of all per-entry rows
    peak_id = 0

    with engine.begin() as connection:
        def flush():
            for model, model_rows in rows.items():
                if model_rows:
                    connection.execute(insert(model.__table__), model_rows)
                    model_rows.clear()

        for library_file in msp_files:
            lib_name = os.path.splitext(os.path.basename(library_file))[0]
            file_size = max(os.path.getsize(library_file), 1)
            n_lib = 0
            for lines_meta, lines_peaks, n_bytes in _iter_entries(library_file):
                # keys are lower case as in MSPReader.create_feature
                entry = {k.lower(): _value_or_none(v) for k, v in _parse_lines(lines_meta).items()}
                if not entry:  # no metadata, e.g. a stray peak line
                    continue
                f_id += 1
                n_lib += 1

                rt_minutes = entry.get('rt_minutes')
                rows[FeatureCombined].append(dict(id=f_id, feature_id=f_id))
                rows[PeakList].append(dict(id=f_id))
                for line in lines_peaks:
                    peak_id += 1
                    mz, intensity, annotation = _parse_peak(line)
                    rows[PeakFeature].append(
                        dict(id=peak_id, peak_list_id=f_id, mz=mz, intensity=intensity, annotation=annotation))
                rows[FeatureMgf].append(dict(id=f_id, feature_id=f_id, combined_feature_id=f_id))
                rows[MsSpec].append(dict(id=f_id, ms_level=2, peaks_id=f_id, feature_mgf_id=f_id))
                rows[FeatureMetaboScape].append(dict(
                    id=f_id,
                    feature_id=f_id,
                    rt_seconds=rt_minutes * 60 if rt_minutes is not None else None,
                    CCS=entry.get('ccs'),
                    mz_meas=entry.get('mz'),
                    adduct_metaboscape=entry.get('ion'),
                    formula_metaboscape=entry.get('formula'),
                    annotation_type=lib_name,  # shown as source of library matches
                    annotation_source=entry.get('comment'),  # abuse annotation source for comment
                    combined_feature_id=f_id,
                ))
                rows[FeatureSirius].append(dict(
                    id=f_id, feature_id=f_id, use_zodiac_scoring_for_best=True, combined_feature_id=f_id))
                rows[CompoundCandidate].append(dict(
                    id=f_id,
                    feature_id=f_id,
                    name_sirius=entry.get('name'),
                    xlogp=entry.get('logp'),
                    inchi=entry.get('inchi'),
                    smiles=entry.get('smiles'),
                    confidence_rank=entry.get('confidence_level'),  # higher rank/level is better
                ))

                if len(rows[FeatureCombined]) >= commit_at_latest_after:
                    flush()
                    if progress is not None:
                        progress(lib_name, n_lib, n_bytes / file_size)
            flush()
            if progress is not None:
                progress(lib_name, n_lib, 1.)
    engine.dispose()
    return f_id

if __name__ == '__main__':
    pass
    # from msIO.feature_managers.db import Library

    # lib = Library(db_file)

