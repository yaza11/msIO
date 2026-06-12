import numpy as np

from msIO import PeakList


def cosine_similarity_forward(ref: PeakList, meas: PeakList, max_dmz_da: float, return_nhits: bool = False) -> tuple[float] | tuple[float, int]:
    """Match peaks of b in a (a is therefore the reference)"""
    if (ref is None) or (meas is None):
        if return_nhits:
            return float('nan'), 0
        return float('nan'),

    norm2 = lambda intensities: sum([i ** 2 for i in intensities]) ** 0.5

    norm_a: float = norm2(ref.intensities)
    norm_b: float = norm2(meas.intensities)

    mzs_ref = np.asarray(ref.mzs, dtype=float)
    ints_ref = np.asarray(ref.intensities, dtype=float)

    running_score: float = 0.
    n_hits: int = 0
    for pb in meas.peaks:
        if np.any(idcs_match := (np.abs(pb.mz - mzs_ref) < max_dmz_da)):
            running_score += pb.intensity * ints_ref[idcs_match].sum()
            n_hits += 1

    score = min(1., running_score / (norm_a * norm_b))
    if return_nhits:
        return score, n_hits
    return score,


def cosine_similarity_backward(ref: PeakList, meas: PeakList, max_dmz_da: float, return_nhits: bool = False) -> tuple[float] | tuple[float, int]:
    return cosine_similarity_forward(meas, ref, max_dmz_da, return_nhits=return_nhits)


def cosine_similarity_sym(a, b, max_dmz_da, return_nhits: bool = False) -> tuple[float] | tuple[float, int]:
    score_fwd, *n_hits_fwd = cosine_similarity_forward(a, b, max_dmz_da, return_nhits=return_nhits)
    score_bwd, *n_hits_bwd = cosine_similarity_forward(b, a, max_dmz_da, return_nhits=return_nhits)
    score = (score_fwd + score_bwd) / 2

    if return_nhits:
        return score, min(n_hits_fwd[0], n_hits_bwd[0])
    return score,


if __name__ == '__main__':
    pl1 = PeakList(mzs=[1, 2, 3], intensities=[1, 2, 3])
    pl2 = PeakList(mzs=[2, 3, 4], intensities=[2, 3, 1])

    score, *n_hits = cosine_similarity_forward(pl1, pl2, 0.5, True)
    score, *n_hits = cosine_similarity_sym(pl1, pl2, 0.5, True)
