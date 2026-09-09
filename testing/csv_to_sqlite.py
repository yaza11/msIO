"""This script shows how to convert a csv table (from mzmine or MetaboScape) into an sqlite file suitable for downstram analysis."""
import os

from msIO import MgfImportManager
from msIO.feature_managers.combined import ProjectImportManager
from msIO.feature_managers.db import FeatureManagerDB
from msIO.feature_managers.gnps import GnpsImportManager
from msIO.feature_managers.metaboscape import MetaboscapeImportManager
from msIO.feature_managers.sirius import SiriusImportManager
from msIO.sql.session import initiate_db

# I/O
path_folder = r'C:\Users\yanni\Nextcloud2\Promotion\msIO\testing\test data'
path_metaboscape_csv = os.path.join(path_folder, "height.csv")
path_mgf = os.path.join(path_folder, "ms1_ms2.mgf")
path_sirius_folder = os.path.join(path_folder, "sirius")
path_gnps_folder = os.path.join(path_folder, "gnps")

db_file = "height.sqlite"

# write the sqlite file
metaboscape_manager = MetaboscapeImportManager(path_metaboscape_csv)
mgf_manager = MgfImportManager(path_mgf, accumulate_spectra=False)
gnps_manager = GnpsImportManager(path_gnps_folder)
sirius_manager = SiriusImportManager(path_folder_export=path_sirius_folder)

project_import_manager = ProjectImportManager(
    metaboscape_manager=metaboscape_manager,
    mgf_manager=mgf_manager,
    gnps_manager=gnps_manager,
    sirius_manager=sirius_manager
)

initiate_db(db_file)
project_import_manager.to_sql(db_file, feature_ids=metaboscape_manager.feature_ids)
