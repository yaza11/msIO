"""
Turn the matches found in the libraries into the results table of the GUI (no GUI code, runs in worker threads).

The results table has one row per match (measured feature, library, library feature) and optionally one dummy row per
measured feature, which only holds the properties of the feature (see dummy_table).
"""

import numpy as np
import pandas as pd
from LipidCalculator import Adduct, CompoundDict

from msIO import PeakList
from msIO.feature_managers.db import FeatureManagerDB, Library, lib_match_result_to_table

# columns of the results table besides the shown ones (see ResultsPanel.COLUMNS)
EXTRA_COLUMNS: list[str] = ['feature_id_lib', 'formula_display', 'default_order', 'is_dummy']
# properties of the measured features, taken from the feature table (see build_feature_table)
FEATURE_COLUMNS: list[str] = ['feature_id_meas', 'mz_meas', 'RT', 'ccs_meas', 'intensity', 'i1_i0']
# library feature id of the dummy rows
DUMMY_LIB_ID: int = -1
# columns of the matches that hold text (empty in dummy rows, the other match columns are missing values)
TEXT_COLUMNS: set[str] = {'name', 'formula', 'formula_display', 'adduct', 'source_library', 'library'}

# distance of the M+1 peak (13C) to the monoisotopic peak of singly charged ions
C13_DISTANCE: float = (CompoundDict('C[13]') - CompoundDict('C')).mass
# M0 and M1 peaks of the MS1 spectrum are searched within this distance of their expected m/z
ISOTOPE_PEAK_TOLERANCE_DA: float = 0.01


def empty_matches_table(shown_columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=list(dict.fromkeys(shown_columns + EXTRA_COLUMNS)))


def isotope_ratio(ms1: PeakList | None, mz: float | None) -> float:
    """
    Intensity ratio of the M+1 and the monoisotopic peak (I1 / I0) in the MS1 spectrum: M0 is the highest peak within
    ISOTOPE_PEAK_TOLERANCE_DA of the precursor m/z, M1 the highest one within the same distance of M0 + 13C.
    NaN if the spectrum has no M0 peak, 0 if it has no M1 peak.
    """
    if ms1 is None or mz is None or len(ms1.mzs) == 0:
        return np.nan
    mzs, intensities = np.asarray(ms1.mzs, dtype=float), np.asarray(ms1.intensities, dtype=float)

    def highest_peak_near(target: float) -> tuple[float, float] | None:
        near = np.abs(mzs - target) <= ISOTOPE_PEAK_TOLERANCE_DA
        if not near.any():
            return None
        idx: int = int(np.flatnonzero(near)[np.argmax(intensities[near])])
        return mzs[idx], intensities[idx]

    if (m0 := highest_peak_near(mz)) is None or m0[1] <= 0:
        return np.nan
    m1 = highest_peak_near(m0[0] + C13_DISTANCE)
    return 0. if m1 is None else m1[1] / m0[1]


def build_feature_table(meas_db: FeatureManagerDB, ms1_spectra: dict[int, PeakList]) -> pd.DataFrame:
    """
    Properties of all measured features shown in the results table: m/z, RT in minutes, CCS, summed intensity and
    I1 / I0. Features without m/z are left out.
    """
    mzs: dict[int, float] = meas_db.mzs
    f_ids: list[int] = sorted(f_id for f_id, mz in mzs.items() if mz is not None)
    table = pd.DataFrame({'feature_id_meas': np.array(f_ids, dtype=int)})
    f_ids_meas: pd.Series = table.feature_id_meas
    table['mz_meas'] = f_ids_meas.map(mzs).astype(float).round(4)
    table['RT'] = (pd.to_numeric(f_ids_meas.map(meas_db.retention_times_in_seconds), errors='coerce') / 60).round(2)
    # CCS of 0 means not measured
    ccs = pd.to_numeric(f_ids_meas.map(meas_db.collisional_cross_sections), errors='coerce')
    table['ccs_meas'] = ccs.where(ccs > 0).round(1)
    table['intensity'] = pd.to_numeric(f_ids_meas.map(meas_db.summed_intensities), errors='coerce').round(2)
    table['i1_i0'] = np.round([isotope_ratio(ms1_spectra.get(f_id), mzs[f_id]) for f_id in f_ids], 3)
    return table


def dummy_table(feature_table: pd.DataFrame, shown_columns: list[str]) -> pd.DataFrame:
    """One row per measured feature without a library match: only the feature properties are filled."""
    table = feature_table.copy()
    for column in empty_matches_table(shown_columns).columns:
        if column not in table.columns:
            table[column] = '' if column in TEXT_COLUMNS else np.nan
    table['feature_id_lib'] = DUMMY_LIB_ID
    table['is_dummy'] = True
    table['default_order'] = np.arange(len(table))
    return table


def combine_with_dummies(matches_table: pd.DataFrame, dummies: pd.DataFrame) -> pd.DataFrame:
    """Matches and dummy rows, the dummy row of a feature comes before its matches in the default order."""
    if matches_table.empty:
        return dummies.copy()
    table = pd.concat([matches_table, dummies], ignore_index=True)
    table['is_dummy'] = table['is_dummy'].astype(bool)
    # stays a whole number with the missing values of the dummy rows (shown as 43, not 43.0)
    table['n_hits_ms2'] = pd.to_numeric(table['n_hits_ms2'], errors='coerce').astype('Int64')
    table = table.sort_values(['feature_id_meas', 'is_dummy', 'default_order'], ascending=[True, False, True],
                              kind='mergesort')
    table['default_order'] = np.arange(len(table))
    return table.reset_index(drop=True)


def limit_matches_per_feature(table: pd.DataFrame, max_matches: int) -> pd.DataFrame:
    """
    Keep at most max_matches rows per measured feature: the dummy row first, then the matches with the highest MS2
    score (matches without score last, ties in the default order).
    """
    if table.empty:
        return table
    ranked = table.assign(_score=pd.to_numeric(table['ms2_score'], errors='coerce').fillna(-np.inf))
    ranked = ranked.sort_values(['feature_id_meas', 'is_dummy', '_score', 'default_order'],
                                ascending=[True, False, False, True], kind='mergesort')
    keep = ranked.groupby('feature_id_meas', sort=False).cumcount() < max_matches
    return table.loc[ranked.index[keep.to_numpy()]].sort_values('default_order', kind='mergesort')


def formula_fits_mz(
        formulas: pd.Series,
        adducts: pd.Series,
        mzs_meas: np.ndarray,
        max_dmz_ppm: float,
        max_dmz_da: float
) -> np.ndarray:
    """
    Whether the m/z calculated from formula and adduct matches the measured m/z within the precursor tolerance
    (stricter of both, as in the search). Formulas or adducts that can't be parsed don't fit.
    """
    def calculate_mz(formula: str, adduct: str) -> float:
        try:
            return Adduct(adduct).mass_to_mz(CompoundDict(formula).mass)
        except Exception:
            return np.nan

    # many matches share formula and adduct
    mzs_calculated: dict[tuple[str, str], float] = {key: calculate_mz(*key) for key in set(zip(formulas, adducts))}
    mzs_calc = np.array([mzs_calculated[key] for key in zip(formulas, adducts)], dtype=float)
    tolerances = np.minimum(max_dmz_da, max_dmz_ppm * 1e-6 * mzs_meas)
    # nan (not parsable) compares as False
    return np.abs(mzs_calc - mzs_meas) <= tolerances


def matches_to_table(
        matches: dict[int, list[dict]],
        feature_table: pd.DataFrame,
        lib_dbs: dict[str, Library],
        max_dmz_ppm: float,
        max_dmz_da: float,
        shown_columns: list[str]
) -> pd.DataFrame:
    """
    Turn the matches of all libraries (output of Library.find_matches, each match with the file name of its library
    under "library") into the table shown in the GUI.

    Adds the properties of the measured features (see build_feature_table), the CCS difference to the library entry
    and the formula to display (in brackets if it doesn't fit the measured m/z).
    Runs in the search worker, so the data is passed in instead of using the current one.
    """
    if not any(matches.values()):  # lib_match_result_to_table can't handle that
        return empty_matches_table(shown_columns)
    table: pd.DataFrame = lib_match_result_to_table(matches)
    table.loc[:, 'ms2_score'] = table['ms2_score'].round(3)
    table.loc[:, 'dmz_mda'] = table['dmz_mda'].round(1).abs()
    table.loc[:, 'dmz_ppm'] = table['dmz_ppm'].round(1).abs()
    table = table.merge(feature_table, on='feature_id_meas', how='left')

    # CCS of 0 means not measured, feature ids of different libraries overlap
    ccs_lib = pd.Series(np.nan, index=table.index)
    for library_name, lib_db in lib_dbs.items():
        rows = table.library == library_name
        ccs_lib[rows] = pd.to_numeric(table.feature_id_lib[rows].map(lib_db.collisional_cross_sections),
                                      errors='coerce')
    ccs_meas = table.ccs_meas
    table['dccs_percent'] = ((ccs_meas - ccs_lib) / ccs_lib * 100).where((ccs_meas > 0) & (ccs_lib > 0)).round(2)

    has_formula = ~table.formula.isin(['', 'None', 'nan'])
    fits_mz = formula_fits_mz(table.formula, table.adduct, table.mz_meas.to_numpy(dtype=float), max_dmz_ppm,
                              max_dmz_da)
    table['formula_display'] = np.where(has_formula, np.where(fits_mz, table.formula, '(' + table.formula + ')'), '')
    table['is_dummy'] = False
    # order of the search result, restored when the sorting is removed (see ResultsPanel)
    table['default_order'] = np.arange(len(table))
    return table
