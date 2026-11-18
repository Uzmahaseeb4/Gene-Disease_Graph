
from pathlib import Path
import pandas as pd

# Inspect data files in the current project folder.
ROOT = Path(__file__).resolve().parent

DATA_EXTENSIONS = {
    ".parquet", ".csv", ".tsv", ".txt",
    ".gz", ".xlsx", ".xls"
}

SKIP_DIRS = {
    ".git", ".venv", "venv", "node_modules",
    "__pycache__", "site-packages"
}


def find_data_files(root):
    for path in root.rglob("*"):
        if not path.is_file():
            continue

        if any(part in SKIP_DIRS for part in path.parts):
            continue

        if any(part.startswith(".") for part in path.relative_to(root).parts):
            continue

        if path.suffix.lower() in DATA_EXTENSIONS:
            yield path


def inspect_file(path):
    print("\n" + "=" * 80)
    print(f"FILE: {path.relative_to(ROOT)}")

    try:
        suffixes = [s.lower() for s in path.suffixes]

        if ".parquet" in suffixes:
            df = pd.read_parquet(path)
        elif path.suffix.lower() in {".xlsx", ".xls"}:
            df = pd.read_excel(path, nrows=5)
        else:
            # Read only a small sample, not the entire CSV/TSV.
            df = pd.read_csv(
                path,
                sep=None,
                engine="python",
                compression="infer",
                nrows=5,
                on_bad_lines="skip",
            )

        print("\nCOLUMNS:")
        print(df.columns.tolist())

        print("\nSAMPLE ROWS:")
        print(df.head(3).to_string(index=False))

    except Exception as exc:
        print(f"Could not inspect this file: {exc}")


def main():
    files = list(find_data_files(ROOT))

    print(f"Project folder: {ROOT}")
    print(f"Data files found: {len(files)}")

    if not files:
        print(
            "No data files found. Check ROOT or the dataset folder location."
        )
        return

    print("\nAVAILABLE FILES:")
    for i, path in enumerate(files, start=1):
        print(f"{i}. {path.relative_to(ROOT)}")

    # Inspect at most 10 files in this first pass.
    for path in files[:10]:
        inspect_file(path)

    if len(files) > 10:
        print(
            f"\n{len(files) - 10} additional files were not inspected. "
            "We can inspect the relevant ones next."
        )


if __name__ == "__main__":
    main()
