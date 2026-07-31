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
    from matchms.similarity import CosineGreedy, NeutralLossesCosine
    from msIO.feature_managers.db import Library, FeatureManagerDB

    path_file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\mzmine\database_all_features.db"
    dbm = FeatureManagerDB(path_file)

    # fetch some example
    # f_id1, f_id2 = 6640, 6641
    # f_id1, f_id2 = 6640, 5835
    # f_id1, f_id2 = 6640, 4558
    f_id1, f_id2 = 5197, 4558
    f1 = dbm.get_ms_spectrum(f_id1, 2)
    f2 = dbm.get_ms_spectrum(f_id2, 2)
    # names = lib.names_metaboscape
    mz1 = dbm.mzs[f_id1]
    mz2 = dbm.mzs[f_id2]

    ax = f1.plot()
    f2.plot(ax=ax, as_mirror=True)

    s1 = peaklist_to_spectrum(f1, precursor_mz=mz1, charge=1)
    s2 = peaklist_to_spectrum(f2, precursor_mz=mz2, charge=1)

    mod = ModifiedCosineGreedy(tolerance=10e-3)
    mod_score = mod.pair(s1, s2)

    cos = CosineGreedy(tolerance=10e-3)
    cos_score = cos.pair(s1, s2)

    neut = NeutralLossesCosine(tolerance=10e-3)
    neut_score = neut.pair(s1, s2)
