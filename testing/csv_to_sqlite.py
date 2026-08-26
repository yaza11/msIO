"""This script shows how to convert a csv table (from mzmine or MetaboScape) into an sqlite file suitable for downstram analysis."""
from msIO.feature_managers.combined import ProjectImportManager
from msIO.feature_managers.db import FeatureManagerDB
from msIO.feature_managers.metaboscape import MetaboscapeImportManager
from msIO.sql.session import initiate_db

# I/O
path_metaboscape_csv = "height.csv"
db_file = "height.sqlite"

# write the sqlite file
mmanager = MetaboscapeImportManager(path_metaboscape_csv)
project_import_manager = ProjectImportManager(metaboscape_manager=mmanager)

initiate_db(db_file)
project_import_manager.to_sql(db_file, feature_ids=project_import_manager.feature_ids)