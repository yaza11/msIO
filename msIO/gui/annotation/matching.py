"""Turn the matches found in the library into the results table of the GUI (no GUI code, runs in a worker thread)."""

import numpy as np
import pandas as pd
from LipidCalculator import Adduct, CompoundDict

from msIO.feature_managers.db import FeatureManagerDB, Library, lib_match_result_to_table

# columns of the results table besides the shown ones (see ResultsPanel.COLUMNS)
EXTRA_COLUMNS: list[str] = ['feature_id_lib', 'formula_display', 'default_order']


def empty_matches_table(shown_columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=shown_columns + EXTRA_COLUMNS)


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
        meas_db: FeatureManagerDB,
        lib_db: Library,
        max_dmz_ppm: float,
        max_dmz_da: float,
        shown_columns: list[str]
) -> pd.DataFrame:
    """
    Turn the output of Library.find_matches into the table shown in the GUI.

    Adds properties of the measured features (m/z, RT in minutes, CCS, summed intensity), the CCS difference
    to the library entry and the formula to display (in brackets if it doesn't fit the measured m/z).
    Runs in the search worker, so the databases are passed in instead of using the current ones.
    """
    if not any(matches.values()):  # lib_match_result_to_table can't handle that
        return empty_matches_table(shown_columns)
    table: pd.DataFrame = lib_match_result_to_table(matches)
    f_ids_meas: pd.Series = table.feature_id_meas
    table.loc[:, 'ms2_score'] = table['ms2_score'].round(3)
    table.loc[:, 'dmz_mda'] = table['dmz_mda'].round(1).abs()
    table.loc[:, 'dmz_ppm'] = table['dmz_ppm'].round(1).abs()
    table['RT'] = (f_ids_meas.map(meas_db.retention_times_in_seconds) / 60).round(2)
    table['intensity'] = pd.to_numeric(f_ids_meas.map(meas_db.summed_intensities), errors='coerce').round(2)

    # CCS of 0 means not measured
    ccs_meas = pd.to_numeric(f_ids_meas.map(meas_db.collisional_cross_sections), errors='coerce')
    ccs_lib = pd.to_numeric(table.feature_id_lib.map(lib_db.collisional_cross_sections), errors='coerce')
    table['ccs_meas'] = ccs_meas.where(ccs_meas > 0).round(1)
    table['dccs_percent'] = ((ccs_meas - ccs_lib) / ccs_lib * 100).where((ccs_meas > 0) & (ccs_lib > 0)).round(2)

    mzs_meas = f_ids_meas.map(meas_db.mzs)
    has_formula = ~table.formula.isin(['', 'None', 'nan'])
    fits_mz = formula_fits_mz(table.formula, table.adduct, mzs_meas.to_numpy(), max_dmz_ppm, max_dmz_da)
    table['formula_display'] = np.where(has_formula, np.where(fits_mz, table.formula, '(' + table.formula + ')'), '')
    table['mz_meas'] = mzs_meas.round(4)
    # order of the search result, restored when the sorting is removed (see ResultsPanel)
    table['default_order'] = np.arange(len(table))
    return table
