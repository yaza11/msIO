import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from pathlib import Path

from msIO.features.base import SqlBaseClass


def get_engine(db_file: str):
    """Return a SQLAlchemy engine for the given SQLite database file."""
    db_path = Path(db_file).resolve()
    return create_engine(f"sqlite+pysqlite:///{db_path}", future=True)


def get_sessionmaker(db_file: str):
    """Return a sessionmaker bound to the SQLite database."""
    engine = get_engine(db_file)
    return sessionmaker(bind=engine, autoflush=False, autocommit=False, future=True)


def initiate_db(path_file):
    # after deleting the database, it has to be reinitialized
    if os.path.exists(path_file):
        os.remove(path_file)

    engine = get_engine(path_file)
    SqlBaseClass.metadata.create_all(engine)


if __name__ == '__main__':
    initiate_db(r"C:\Users\Yannick Zander\Nextcloud2\Promotion\msIO\msIO\feature_managers\database.db")
