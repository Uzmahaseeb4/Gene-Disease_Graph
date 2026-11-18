
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
IMPORT_DIR = ROOT / "data" / "graph" / "neo4j_import"

FILES = [
    "genes.csv",
    "diseases.csv",
    "associations.csv",
]

def inspect_file(filename):
    path = IMPORT_DIR / filename

    print("\n" + "=" * 70)
    print(f"FILE: {path}")

    if not path.exists():
        print("File not found!")
        return

    print(f"File size: {path.stat().st_size / (1024 * 1024):.2f} MB")

    try:
        df = pd.read_csv(path, nrows=5)

        print("\nColumns:")
        print(df.columns.tolist())

        print("\nSample records:")
        print(df.head(5).to_string(index=False))

        print("\nMissing values in sample:")
        print(df.isna().sum().to_string())

    except Exception as exc:
        print(f"Could not inspect file: {exc}")

def main():
    for filename in FILES:
        inspect_file(filename)

if __name__ == "__main__":
    main()

