"""
Series of measured features (e.g. homologues differing by CH2): fit a line through features selected in the scatter
plot and extrapolate it (no GUI code).

Members of a series are spaced by a constant m/z step, member n of the series has the m/z mz_0 + n * step and the
RT rt_0 + n * rt_step. The step is estimated from the selected features, which don't have to be consecutive members.
"""

from dataclasses import dataclass, field

import numpy as np
from LipidCalculator import CompoundDict

# m/z differences of features closer than this are not used to estimate the step (e.g. duplicates)
MIN_STEP_DA: float = 0.5
# the m/z differences of the selected features have to be multiples of the step within this fraction of the step
MAX_STEP_DEVIATION: float = 0.1
# the step is at most this many times smaller than the smallest m/z difference of the selected features
MAX_STEP_DIVISOR: int = 6
# common steps of series, the closest one within KNOWN_STEP_TOLERANCE_DA is shown with the fitted step
KNOWN_STEPS: dict[str, float] = {name: CompoundDict(name).mass for name in (
    'H2', 'CH2', 'C2H4', 'O', 'H2O', 'C2H2', 'C5H8', 'C3H6', 'C6H10O5'
)}
KNOWN_STEP_TOLERANCE_DA: float = 0.005


@dataclass
class SeriesFit:
    """
    Line through the selected features of a series.

    members: member index n of each selected feature (same order as the features), the lowest is 0
    n_min, n_max: members shown (selected ones and the extrapolated ones)
    max_dmz_da, max_drt: largest deviation of the selected features from the line
    """
    step: float
    mz_0: float
    rt_step: float
    rt_0: float
    members: np.ndarray
    n_min: int
    n_max: int
    max_dmz_da: float
    max_drt: float
    step_name: str | None = None
    # predicted members (n, rt, mz) without a selected feature
    predicted: list[tuple[int, float, float]] = field(default_factory=list)

    def mz(self, n) -> np.ndarray:
        return self.mz_0 + np.asarray(n) * self.step

    def rt(self, n) -> np.ndarray:
        return self.rt_0 + np.asarray(n) * self.rt_step

    def describe(self) -> str:
        name: str = f" ({self.step_name})" if self.step_name else ""
        return f"Δm/z {self.step:.4f}{name}, ΔRT {self.rt_step:+.2f} min"


def estimate_step(mzs: np.ndarray) -> float | None:
    """
    m/z step of a series through mzs: the smallest m/z difference (above MIN_STEP_DA) or the largest fraction of it of
    which all differences are multiples, e.g. features at 0, 28 and 42 Da give 14 Da. The step is at least
    MIN_STEP_DA. None if no step fits.
    """
    differences: np.ndarray = np.diff(np.sort(mzs))
    differences = differences[differences >= MIN_STEP_DA]
    if differences.size == 0:
        return None
    smallest: float = float(differences.min())
    offsets: np.ndarray = np.sort(mzs) - np.min(mzs)
    for divisor in range(1, MAX_STEP_DIVISOR + 1):
        if (step := smallest / divisor) < MIN_STEP_DA:
            break
        multiples: np.ndarray = offsets / step
        if np.all(np.abs(multiples - np.round(multiples)) <= MAX_STEP_DEVIATION):
            return step
    return None


def name_step(step: float) -> str | None:
    """Name of the known step closest to step (within KNOWN_STEP_TOLERANCE_DA)."""
    name, mass = min(KNOWN_STEPS.items(), key=lambda item: abs(item[1] - step))
    return name if abs(mass - step) <= KNOWN_STEP_TOLERANCE_DA else None


def fit_series(rts: np.ndarray, mzs: np.ndarray, n_extrapolate: int) -> SeriesFit | None:
    """
    Fit a series through the features at (rts, mzs): estimate the step, number the members and fit m/z and RT
    linearly in the member number. The line covers the selected features and n_extrapolate members beyond the first
    and the last one. None if the features don't define a series (less than two different m/z or no common step).
    """
    rts, mzs = np.asarray(rts, dtype=float), np.asarray(mzs, dtype=float)
    if (step := estimate_step(mzs)) is None:
        return None
    members: np.ndarray = np.round((mzs - mzs.min()) / step).astype(int)
    # least squares with the member numbers: mz = mz_0 + n * step, rt = rt_0 + n * rt_step
    step, mz_0 = np.polyfit(members, mzs, 1)
    rt_step, rt_0 = np.polyfit(members, rts, 1)
    fit = SeriesFit(
        step=float(step), mz_0=float(mz_0), rt_step=float(rt_step), rt_0=float(rt_0), members=members,
        n_min=int(members.min()) - n_extrapolate, n_max=int(members.max()) + n_extrapolate,
        max_dmz_da=float(np.max(np.abs(mz_0 + members * step - mzs))),
        max_drt=float(np.max(np.abs(rt_0 + members * rt_step - rts))),
        step_name=name_step(float(step)),
    )
    occupied: set[int] = set(members.tolist())
    fit.predicted = [(n, float(fit.rt(n)), float(fit.mz(n))) for n in range(fit.n_min, fit.n_max + 1)
                     if n not in occupied and fit.mz(n) > 0]
    return fit


def find_candidates(
        fit: SeriesFit,
        f_ids: np.ndarray,
        rts: np.ndarray,
        mzs: np.ndarray,
        mz_tolerances: np.ndarray,
        rt_window: float
) -> dict[int, list[int]]:
    """
    Measured features at the predicted members of a series: {member n: [feature ids]}, matched with the m/z tolerance
    of each feature and the RT window around the predicted RT.
    """
    candidates: dict[int, list[int]] = {}
    for n, rt, mz in fit.predicted:
        near = (np.abs(mzs - mz) <= mz_tolerances) & (np.abs(rts - rt) <= rt_window)
        if near.any():
            candidates[n] = [int(f_id) for f_id in f_ids[near]]
    return candidates
