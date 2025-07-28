import os

from msIO.feature_managers.main import SIRIUS_FILE_NAMES, get_sirius_file_for_tag


class ExportManager:
    metaboscape_file: str = None
    gnps_file: str = None

    def __init__(
            self,
            path_metaboscape_folder: str,
            path_gnps_folder: str,
            path_sirius_export_folder: str,
            sirius_export_tag: str = None
    ) -> None:
        self.path_metaboscape_folder = path_metaboscape_folder
        self.path_gnps_folder = path_gnps_folder
        self.path_sirius_export_folder = path_sirius_export_folder
        self.sirius_export_tag = sirius_export_tag

    def _find_metaboscape_files(self):
        for file in os.listdir(self.path_metaboscape_folder):
            if file.split('.', 1)[1] == 'csv':
                self.metaboscape_file = file
                break

    def _find_gnps_file(self):
        for file in os.path.join(self.path_gnps_folder, 'gnps_molecular_network_graphml'):
            if file.endswith('.graphml'):
                self.gnps_file = file
                break


if __name__ == '__main__':
    em = ExportManager(
        path_metaboscape_folder=r'\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas\U1545B_U1549B\MetabSscape'
    )
