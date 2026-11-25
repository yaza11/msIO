import numpy as np
from matchms import Spectrum
from matchms.similarity import ModifiedCosineGreedy

from msIO import PeakList


def peaklist_to_spectrum(
        peaklist: PeakList | None,
        precursor_mz: float = None,
        charge: int = None,
) -> Spectrum | None:
    if peaklist is None:
        return None

    mzs = np.array(peaklist.mzs)
    intensities = np.array(peaklist.intensities)

    # Filter out zero/negative intensities and m/z <= 0
    valid = (mzs > 0) & (intensities > 0)
    mzs = mzs[valid]
    intensities = intensities[valid]

    # Sort by m/z (required by matchms)
    idx = np.argsort(mzs)
    mzs = mzs[idx]
    intensities = intensities[idx]

    # Build metadata
    metadata = {}
    if precursor_mz is not None:
        metadata["precursor_mz"] = float(precursor_mz)
    if charge is not None:
        metadata["charge"] = int(charge)

    return Spectrum(mz=mzs, intensities=intensities, metadata=metadata)



if __name__ == '__main__':
    pass

