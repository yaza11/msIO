from sqlalchemy import create_engine, MetaData, text
import subprocess
import sys
import os


def dump_schema(db_url: str, output_file: str):
    from sqlalchemy import create_engine, MetaData, text

    engine = create_engine(db_url, future=True)
    metadata = MetaData()
    metadata.reflect(bind=engine)

    create_statements = []
    with engine.connect() as conn:
        for table_name in metadata.tables.keys():
            if engine.dialect.name == "mysql":
                result = conn.execute(text(f"SHOW CREATE TABLE `{table_name}`"))
                row = result.fetchone()
                if row:
                    create_statements.append(f"{row[1]};\n")
            elif engine.dialect.name == "sqlite":
                result = conn.execute(
                    text("SELECT sql FROM sqlite_master WHERE type='table' AND name=:name"),
                    {"name": table_name}
                )
                row = result.fetchone()
                if row and row[0]:
                    create_statements.append(f"{row[0]};\n")
            else:
                raise NotImplementedError(f"Unsupported DB dialect: {engine.dialect.name}")

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n\n".join(create_statements))

    print(f"[OK] {engine.dialect.name} schema dumped to: {output_file}")


def generate_er_diagram(db_url: str, diagram_file: str):
    """
    Generate an ER diagram from the database using ERAlchemy.
    """
    try:
        subprocess.run(
            ["eralchemy", "-i", db_url, "-o", diagram_file],
            check=True
        )
        print(f"[OK] Diagram generated at: {diagram_file}")
    except FileNotFoundError:
        print("[ERROR] ERAlchemy not installed. Install with:")
        print("    pip install eralchemy")
    except subprocess.CalledProcessError as e:
        print("[ERROR] Failed to generate diagram:", e)


if __name__ == "__main__":
    # ---- USER CONFIG ----
    # Example DB: mysql+pymysql://user:pass@localhost:3306/dbname
    DB_URL = r"sqlite+pysqlite:///C:\Users\Yannick Zander\Nextcloud2\Promotion\msIO\msIO\feature_managers\database.db"
    SQL_FILE = r"mysql_schema.sql"
    DIAGRAM_FILE = r"mysql_schema.pdf"
    # ---------------------

    dump_schema(DB_URL, SQL_FILE)
    # generate_er_diagram(DB_URL, DIAGRAM_FILE)
