import numpy as np

from msIO import PeakList


def cosine_similarity_forward(ref: PeakList, meas: PeakList, max_dmz_da: float) -> float:
    """Match peaks of b in a (a is therefore the reference)"""
    if (ref is None) or (meas is None):
        return 0.

    norm2 = lambda intensities: sum([i ** 2 for i in intensities])

    norm_a = norm2(ref.intensities)
    norm_b = norm2(meas.intensities)

    mzs_ref = np.asarray(ref.mzs, dtype=float)
    ints_ref = np.asarray(ref.intensities, dtype=float)

    running_score = 0
    for pb in meas.peaks:
        if np.any(idcs_match := (np.abs(pb.mz - mzs_ref) < max_dmz_da)):
            running_score += pb.intensity * ints_ref[idcs_match].sum()

    return running_score / (norm_a * norm_b)


def cosine_similarity_backward(ref: PeakList, meas: PeakList, max_dmz_da: float) -> float:
    return cosine_similarity_forward(meas, ref, max_dmz_da)


def cosine_similarity_sim(a, b, max_dmz_da) -> float:
    score = (
            cosine_similarity_forward(a, b, max_dmz_da) / 2
            + cosine_similarity_forward(b, a, max_dmz_da) / 2
    )
    return score
