from msIO.feature_managers.db import FeatureManagerDB, Library

path_file_meas = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\mzmine\database_all_features.db"
path_file_lib = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\compounds\sql\library_complete.sql"


meas = FeatureManagerDB(path_file_meas)
lib = Library(path_file_lib)

