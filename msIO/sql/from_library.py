import os

from sqlalchemy import update
from tqdm import tqdm

from msIO import MSPReader
from msIO.features.combined import FeatureCombined
from msIO.features.metaboscape import FeatureMetaboScape
from msIO.features.mgf import FeatureMgf
from msIO.features.sirius import FeatureSirius
from msIO.sql.session import initiate_db, get_sessionmaker


def write_lib_from_msp_files(db_file: str, msp_files: list[str], commit_at_latest_after=10_000) -> None:
    initiate_db(db_file)

    Session = get_sessionmaker(db_file)
    f_id = 1
    with Session() as session:
        uncommited_features = 0
        for library_file in library_files:
            low_memory = os.path.getsize(library_file) > 1024 ** 3  # 1 GB
            # low_memory = True
            msp_manager = MSPReader(library_file, splitter_peaks_list=None, low_memory=low_memory)
            lib_name = library_file.split('\\')[-1].split('.')[0]

            for i in tqdm(
                    range(msp_manager.n_features),
                    desc='adding features to DB',
                    total=msp_manager.n_features,
                    smoothing=1 / 50
            ):
                if low_memory:
                    msp_manager.read_next()
                    if msp_manager.df_features.empty:
                        continue
                    idx = 0
                else:
                    idx = msp_manager.df_features.index[i]
                f = msp_manager.create_feature(idx, f_id)
                f.metaboscape.annotation_type = lib_name
                session.add(f)
                f_id += 1
                uncommited_features += 1
                if uncommited_features >= COMMIT_AT_LEAST_EVERY:
                    session.commit()
                    uncommited_features = 0
        session.commit()


if __name__ == '__main__':
    # from msIO.feature_managers.db import Library

    # lib = Library(db_file)

    db_file = r"C:\Users\yanni\Downloads\library_high_confidence.sqlite"
    # db_file =  r"C:\Users\Yannick Zander\Downloads\library_high_confidence.sqlite"
    folder_julius = r'\\hlabstorage.dmz.marum.de\scratch\Yannick\compounds\julius\fragments'

    library_files = [
                        # r"C:\Users\Yannick Zander\Downloads\Archlips_High_confidence_spectral_library.msp",
                        # r"C:\Users\Yannick Zander\Downloads\MSMS-Public_experimentspectra-pos-VS19.msp",
                        # r"C:\Users\Yannick Zander\Downloads\Archlips_M+1_Full_spectral_library.msp",
                        # r"C:\Users\yanni\Downloads\Archlips_M+1_Full_spectral_library.msp",
                        r"C:\Users\yanni\Downloads\Archlips_High_confidence_spectral_library.msp",
                        r"C:\Users\yanni\Downloads\MSMS-Public_experimentspectra-pos-VS19.msp"
                        # r"C:\Users\yanni\Downloads\Archlips_Full_spectral_library.msp"
                    ] + [
                        os.path.join(folder_julius, f)
                        for f in os.listdir(folder_julius)
                        if (not f.startswith('!') and not ('test_as_H+' in f))
                    ]

