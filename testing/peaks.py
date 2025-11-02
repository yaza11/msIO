from msIO.list_of_ions.base import PeakFeature, PeakList

peak = PeakFeature(mz=100, intensity=10)

pl = PeakList(mzs=[1, 2, 3], intensities=[3, 4, 5])

pl.peaks.append(peak)
