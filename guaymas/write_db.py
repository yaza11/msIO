"""this now lives in the Guaymas project"""
from msIO import MgfImportManager
from guaymas.paths import path_metaboscape_csv, path_mgf_sirius, \
    path_sirius_folder, db_file, path_gnps_folder
from msIO.feature_managers.combined import ProjectImportManager

from msIO.feature_managers.gnps import GnpsImportManager
from msIO.feature_managers.metaboscape import MetaboscapeImportManager
from msIO.feature_managers.sirius import SiriusImportManager

metaboscape = MetaboscapeImportManager(path_metaboscape_csv)
mgf = MgfImportManager(path_mgf_sirius)
# mgf = None
gnps = GnpsImportManager(path_gnps_folder=path_gnps_folder)
sr = SiriusImportManager(path_folder_export=path_sirius_folder, export_tag='all')

project_import_manager = ProjectImportManager(mgf_manager=mgf,
                                              gnps_manager=gnps,
                                              sirius_manager=sr,
                                              metaboscape_manager=metaboscape)

from msIO.sql.session import initiate_db

initiate_db(db_file)
project_import_manager.to_sql(db_file, feature_ids=project_import_manager.feature_ids[:10])
# project_import_manager.to_sql(db_file)


f = project_import_manager.active_managers['metaboscape'].get_feature(1)

