import os

from tqdm import tqdm

from msIO import MSPReader
from msIO.sql.session import initiate_db, get_sessionmaker


COMMIT_AT_LEAST_EVERY = 1_000

db_file =  r"C:\Users\yanni\Downloads\library_complete.sql"
folder_julius = r'\\hlabstorage.dmz.marum.de\scratch\Yannick\compounds\julius\fragments'

library_files = [
    # r"C:\Users\yanni\Downloads\Archlips_M+1_Full_spectral_library.msp",
    r"C:\Users\yanni\Downloads\Archlips_High_confidence_spectral_library.msp",
    r"C:\Users\yanni\Downloads\MSMS-Public_experimentspectra-pos-VS19.msp"
    # r"C:\Users\yanni\Downloads\Archlips_Full_spectral_library.msp"
] + [
    os.path.join(folder_julius, f)
    for f in os.listdir(folder_julius)
    if (not f.startswith('!') and not ('test_as_H+' in f))
]

# initiate_db(db_file)

Session = get_sessionmaker(db_file)
with Session() as session:
    uncommited_features = 0
    for library_file in library_files:
        msp_manager = MSPReader(library_file, splitter_peaks_list=None)
        feature_ids = msp_manager.df_features.index

        for f_id in tqdm(
                feature_ids,
                desc='adding features to DB',
                total=len(feature_ids),
                smoothing=1/50
        ):
            f = msp_manager.get_feature(f_id)
            session.add(f)
            uncommited_features += 1
            if uncommited_features >= COMMIT_AT_LEAST_EVERY:
                session.commit()
                uncommited_features = 0
    session.commit()

if __name__ == '__main__':
    from msIO.feature_managers.db import FeatureManagerDB

    db_m = FeatureManagerDB(db_file)
