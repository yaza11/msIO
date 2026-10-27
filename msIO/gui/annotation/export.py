"""Write annotated measured features to a new database that can be used as library (no GUI code)."""

import os

from msIO import PeakList
from msIO.feature_managers.db import FeatureManagerDB, Library
from msIO.features.combined import FeatureCombined
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.mgf import FeatureMgf, MsSpec
from msIO.features.sirius import CompoundCandidate, FeatureSirius
from msIO.gui.annotation.common import AnnotationKey
from msIO.sql.session import get_sessionmaker, initiate_db

# kinds of annotations that describe a compound (in-source fragments only refer to another feature)
EXPORTED_KINDS: tuple[str, ...] = ("library", "manual")


def _copy_spectrum(spectrum: PeakList | None) -> PeakList | None:
    """New peak list for the library database (keeps the peak annotations)."""
    if spectrum is None or len(spectrum.mzs) == 0:
        return None
    return PeakList(mzs=list(spectrum.mzs), intensities=list(spectrum.intensities),
                    annotations=list(spectrum.annotations))


def _complete_from_library(kind: str, info: dict, lib_dbs: dict[str, Library]) -> dict:
    """Formula and structure of library annotations (not stored with older annotations) from the loaded library."""
    values: dict = dict(name=info.get('name'), formula=info.get('formula'), adduct=info.get('adduct'),
                        smiles=info.get('smiles'), inchi=info.get('inchi'))
    library_name = info.get('library_name')
    if library_name is None and info.get('library') is not None:
        library_name = os.path.basename(info['library'].path_file)
    if kind == "library" and (lib_db := lib_dbs.get(library_name)) is not None:
        f_id_lib: int = int(info['reference'])
        values['formula'] = values['formula'] or lib_db.formula_metaboscape.get(f_id_lib)
        values['adduct'] = values['adduct'] or lib_db.adduct.get(f_id_lib)
        values['smiles'] = values['smiles'] or lib_db.smiles.get(f_id_lib)
        values['inchi'] = values['inchi'] or lib_db.inchis.get(f_id_lib)
    return values


def write_annotations_as_library(
        path: str,
        meas_db: FeatureManagerDB,
        annotations: dict[AnnotationKey, dict],
        lib_dbs: dict[str, Library],
        source_name: str
) -> int:
    """
    Write one library entry per library match or manual annotation of a measured feature to a new database at path
    (an existing file is replaced). An entry holds the compound of the annotation (name, formula, adduct, structure)
    and the measured feature (m/z, RT, CCS, MS1 and MS2 spectra with their peak annotations), like the entries
    created from msp files (see MSPReader.create_feature). source_name is shown as source of the entries.

    Returns the number of entries.
    """
    entries = sorted(((key, info) for key, info in annotations.items() if key[0] in EXPORTED_KINDS),
                     key=lambda item: (item[0][1], item[0][0], item[0][2]))
    f_ids: list[int] = sorted({key[1] for key, _ in entries})
    ms1_spectra: dict[int, PeakList] = meas_db.get_ms_spectra(f_ids, level=1)
    ms2_spectra: dict[int, PeakList] = meas_db.get_ms_spectra(f_ids, level=2)
    meas_name: str = os.path.basename(meas_db.path_file)

    initiate_db(path)
    with get_sessionmaker(path)() as session:
        for f_id_new, ((kind, f_id_meas, _), info) in enumerate(entries, start=1):
            values: dict = _complete_from_library(kind, info, lib_dbs)
            mz: float | None = meas_db.mzs.get(f_id_meas)
            rt: float | None = meas_db.retention_times_in_seconds.get(f_id_meas)
            ccs: float | None = meas_db.collisional_cross_sections.get(f_id_meas)
            ms_specs = [MsSpec(peaks=peaks, ms_level=level, mz=mz, rt_seconds=rt, ion=values['adduct'])
                        for level, peaks in ((1, _copy_spectrum(ms1_spectra.get(f_id_meas))),
                                             (2, _copy_spectrum(ms2_spectra.get(f_id_meas))))
                        if peaks is not None]
            origin: str = f"{kind} annotation of feature {f_id_meas} in {meas_name}"
            if kind == "library":
                origin += f" (library feature {info.get('reference')} of {info.get('library_name') or 'library'})"
            session.add(FeatureCombined(
                feature_id=f_id_new,
                metaboscape=FeatureMetaboScape(
                    feature_id=f_id_new,
                    rt_seconds=rt,
                    CCS=ccs if ccs is not None and ccs > 0 else None,
                    mz_meas=mz,
                    adduct_metaboscape=values['adduct'],
                    formula_metaboscape=values['formula'],
                    annotation_type=source_name,  # shown as source of library matches
                    annotation_source=origin,  # comment, as for msp libraries
                ),
                mgf=FeatureMgf(feature_id=f_id_new, mz=mz, rt_seconds=rt, ion=values['adduct'], ms_specs=ms_specs),
                sirius=FeatureSirius(feature_id=f_id_new, compound_candidates=[CompoundCandidate(
                    feature_id=f_id_new,
                    name_sirius=values['name'],
                    formula_sirius=values['formula'],
                    adduct_sirius=values['adduct'],
                    smiles=values['smiles'],
                    inchi=values['inchi'],
                    confidence_score=info.get('ms2_score'),
                )]),
            ))
        session.commit()
    return len(entries)
