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
    found = features_near([(rt, mz) for _, rt, mz in fit.predicted], f_ids, rts, mzs, mz_tolerances, rt_window)
    return {fit.predicted[i][0]: f_ids_i for i, f_ids_i in found.items()}


def features_near(
        points: list[tuple[float, float]],
        f_ids: np.ndarray,
        rts: np.ndarray,
        mzs: np.ndarray,
        mz_tolerances: np.ndarray,
        rt_window: float
) -> dict[int, list[int]]:
    """Measured features at points (rt, mz): {index of the point: [feature ids]}, points without features left out."""
    found: dict[int, list[int]] = {}
    for i, (rt, mz) in enumerate(points):
        near = (np.abs(mzs - mz) <= mz_tolerances) & (np.abs(rts - rt) <= rt_window)
        if near.any():
            found[i] = [int(f_id) for f_id in f_ids[near]]
    return found


@dataclass
class SeriesMesh:
    """
    Mesh spanned by two series: through each member of one series runs a line parallel to the other series, so all
    members (selected, interpolated and extrapolated) are knots of the mesh.

    lines_a: lines parallel to series a (one through each member of series b), lines_b the other way round, each line
    as ((rt, mz), (rt, mz)). knots: intersections of the lines (rt, mz).
    """
    lines_a: list[tuple[tuple[float, float], tuple[float, float]]]
    lines_b: list[tuple[tuple[float, float], tuple[float, float]]]
    knots: list[tuple[float, float]]


# the members of two series closest to their intersection are the same (shared) member if they are this close
SHARED_MEMBER_TOLERANCE_DA: float = 0.01
SHARED_MEMBER_TOLERANCE_MIN: float = 0.2


def mesh_of_series(fit_a: SeriesFit, fit_b: SeriesFit) -> SeriesMesh | None:
    """
    Mesh of two series (see SeriesMesh), None if they are parallel. Each line covers the members of the other series
    and the intersection of both series.

    If the members of both series closest to their intersection are the same (within SHARED_MEMBER_TOLERANCE_DA and
    _MIN, e.g. a feature selected in both series), the mesh starts at the mean of both: the intersection of the fitted
    lines can be far off in m/z if they cross at a small angle (e.g. a series without RT shift). Otherwise the mesh
    starts at the intersection and contains both series as lines.
    """
    origin_a = np.array([fit_a.rt_0, fit_a.mz_0])
    origin_b = np.array([fit_b.rt_0, fit_b.mz_0])
    v_a = np.array([fit_a.rt_step, fit_a.step])
    v_b = np.array([fit_b.rt_step, fit_b.step])
    # intersection: origin_a + s * v_a = origin_b + t * v_b, RT and m/z have different scales, so compare the
    # determinant with the lengths of both vectors
    matrix = np.column_stack([v_a, -v_b])
    if abs(np.linalg.det(matrix)) < 1e-9 * np.linalg.norm(v_a) * np.linalg.norm(v_b):
        return None
    s, t = np.linalg.solve(matrix, origin_b - origin_a)
    members_a = np.arange(fit_a.n_min, fit_a.n_max + 1)
    members_b = np.arange(fit_b.n_min, fit_b.n_max + 1)
    member_a, member_b = origin_a + round(s) * v_a, origin_b + round(t) * v_b
    if abs(member_a[1] - member_b[1]) <= SHARED_MEMBER_TOLERANCE_DA and             abs(member_a[0] - member_b[0]) <= SHARED_MEMBER_TOLERANCE_MIN:
        # positions relative to the shared member in steps
        origin = (member_a + member_b) / 2
        coords_a = members_a - round(s)
        coords_b = members_b - round(t)
    else:
        # positions relative to the intersection in steps (member n of a is at n - s, of b at n - t), 0 is the other
        # series itself, so the members of both series are knots
        origin = origin_a + s * v_a
        coords_a = np.sort(np.append(members_a - s, 0.))
        coords_b = np.sort(np.append(members_b - t, 0.))
    range_a = coords_a.min(), coords_a.max()
    range_b = coords_b.min(), coords_b.max()

    def point(a: float, b: float) -> tuple[float, float]:
        rt, mz = origin + a * v_a + b * v_b
        return float(rt), float(mz)

    return SeriesMesh(
        lines_a=[(point(range_a[0], b), point(range_a[1], b)) for b in coords_b],
        lines_b=[(point(a, range_b[0]), point(a, range_b[1])) for a in coords_a],
        knots=[point(a, b) for a in coords_a for b in coords_b],
    )


@dataclass
class ParallelSeries:
    """
    Series b runs parallel to series a, e.g. the same compounds as another adduct: their coeluting features all differ
    by the same m/z. dmz, drt: mean differences (b - a) of the coeluting features, n_pairs: their number,
    name: adduct exchange that fits dmz (e.g. "[M+H]+ → [M+NH4]+").
    """
    index_a: int
    index_b: int
    dmz: float
    drt: float
    n_pairs: int
    name: str | None = None

    def describe(self) -> str:
        name: str = f" ({self.name})" if self.name else ""
        return (f"Series {self.index_b + 1} parallel to series {self.index_a + 1}: Δm/z {self.dmz:+.4f}{name}, "
                f"ΔRT {self.drt:+.2f} min ({self.n_pairs} coeluting pair{'s' if self.n_pairs > 1 else ''})")


def name_adduct_exchange(dmz: float, adducts: list[str], tolerance: float) -> str | None:
    """Exchange of singly charged adducts of one molecule (e.g. [M+H]+ → [M+NH4]+) whose m/z differ by dmz."""
    from LipidCalculator import Adduct
    best: tuple[float, str] | None = None
    for name_a in adducts:
        for name_b in adducts:
            try:
                adduct_a, adduct_b = Adduct(name_a), Adduct(name_b)
            except Exception:
                continue
            if name_a == name_b or adduct_a.charge != adduct_b.charge or abs(adduct_a.charge) != 1 \
                    or adduct_a.multiplicity != 1 or adduct_b.multiplicity != 1:
                continue
            deviation = abs(adduct_b.mass_to_mz(1000.) - adduct_a.mass_to_mz(1000.) - dmz)
            if deviation <= tolerance and (best is None or deviation < best[0]):
                best = (deviation, f"{name_a} → {name_b}")
    return None if best is None else best[1]


def find_parallel(
        fit_a: SeriesFit,
        points_a: np.ndarray,
        fit_b: SeriesFit,
        points_b: np.ndarray,
        rt_window: float,
        mz_tolerance: float
) -> tuple[float, float, int] | None:
    """
    Whether series b runs parallel to series a: same m/z step and the same m/z difference (within twice the m/z
    tolerance) for all pairs of coeluting features (within the RT window), points are (rt, mz) of the selected
    features. Returns the mean m/z and RT difference (b - a) and the number of pairs, None if not parallel (or no
    features coelute).
    """
    if abs(fit_a.step - fit_b.step) > 2 * mz_tolerance:
        return None
    pairs = [(rt_b - rt_a, mz_b - mz_a) for rt_a, mz_a in points_a for rt_b, mz_b in points_b
             if abs(rt_b - rt_a) <= rt_window and abs(mz_b - mz_a) > mz_tolerance]
    if not pairs:
        return None
    drts, dmzs = np.array(pairs).T
    if dmzs.max() - dmzs.min() > 2 * mz_tolerance:
        return None
    return float(dmzs.mean()), float(drts.mean()), len(pairs)


def shift_mesh(mesh: SeriesMesh, drt: float, dmz: float) -> SeriesMesh:
    """The mesh moved by drt and dmz (e.g. to the other adduct of parallel series)."""
    def shift(point: tuple[float, float]) -> tuple[float, float]:
        return point[0] + drt, point[1] + dmz
    return SeriesMesh(lines_a=[(shift(p), shift(q)) for p, q in mesh.lines_a],
                      lines_b=[(shift(p), shift(q)) for p, q in mesh.lines_b],
                      knots=[shift(p) for p in mesh.knots])


@dataclass
class MeshOfSeries:
    """
    Mesh spanned by the crossing series index_a and index_b (meshes[0]) and its copies (meshes[1:]) for the series
    parallel to one of them (copy_of[k]: index of the parallel series of meshes[k + 1]), e.g. other adducts.
    """
    index_a: int
    index_b: int
    meshes: list[SeriesMesh]
    copy_of: list[int] = field(default_factory=list)

    @property
    def knots(self) -> list[tuple[float, float]]:
        return [knot for mesh in self.meshes for knot in mesh.knots]

    def describe(self) -> str:
        text: str = f"Mesh of series {self.index_a + 1} and {self.index_b + 1}"
        if self.copy_of:
            text += f" (copies for series {', '.join(str(k + 1) for k in self.copy_of)})"
        return text


def mesh_with_parallels(fits: list[SeriesFit | None], parallels: list[ParallelSeries]) -> MeshOfSeries | None:
    """
    Mesh of the two most recent crossing series (not parallel to each other), copied by the m/z and RT offset to each
    series parallel to one of them (e.g. the same series as another adduct), so a parallel series needs no crossing
    series of its own. None if no two fitted series cross.
    """
    parallel_pairs: set[tuple[int, int]] = {(p.index_a, p.index_b) for p in parallels}
    for b in range(len(fits) - 1, 0, -1):
        for a in range(b - 1, -1, -1):
            if fits[a] is None or fits[b] is None or (a, b) in parallel_pairs:
                continue
            if (base := mesh_of_series(fits[a], fits[b])) is None:
                continue
            result = MeshOfSeries(index_a=a, index_b=b, meshes=[base])
            for p in parallels:
                # offset from the base series to the parallel one
                if p.index_a in (a, b) and p.index_b not in (a, b) and p.index_b not in result.copy_of:
                    copy, drt, dmz = p.index_b, p.drt, p.dmz
                elif p.index_b in (a, b) and p.index_a not in (a, b) and p.index_a not in result.copy_of:
                    copy, drt, dmz = p.index_a, -p.drt, -p.dmz
                else:
                    continue
                result.meshes.append(shift_mesh(base, drt, dmz))
                result.copy_of.append(copy)
            return result
    return None
