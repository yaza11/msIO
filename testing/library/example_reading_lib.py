from LipidCalculator import CompoundDict
from LipidCalculator.adduct.parser import Adduct

from msIO.feature_managers.db import Library

path_file_lib = r"..\test data\lib\archlipids_high_conf.sqlite"
lib = Library(path_file_lib)

f_id1 = lib.find_by_name('OH-GTGT-0a')[0]
f_id2 = lib.find_by_name('OH-GTGT-0b')[0]

lib.plot_compound_overview(f_id1)
lib.plot_compound_overview(f_id2)

