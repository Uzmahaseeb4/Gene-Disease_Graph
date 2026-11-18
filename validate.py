
from pathlib import Path
import hashlib
import pandas as pd

CLEAN_DIR = Path("data/clean")
RAW_DIR = Path("data/raw")

FILES = [
    "genes.tsv",
    "diseases.tsv",
    "disease_hierarchy.tsv",
    "disease_therapeutic_area.tsv",
    "gene_disease_associations.tsv",
    "gwas_variants.tsv",
    "gwas_studies.tsv",
    "gwas_variant_study.tsv",
    "gwas_study_disease.tsv",
    "gwas_variant_gene.tsv",
]

def normalize(name):
    """Make column names easier to compare."""
    return "".join(ch.lower() for ch in str(name) if ch.isalnum())

def find_column(df, candidates):
    """Find a column using possible alternative names."""
    lookup = {normalize(c): c for c in df.columns}
    for candidate in candidates:
        found = lookup.get(normalize(candidate))
        if found is not None:
            return found
    return None

def read_tsv(filename):
    path = CLEAN_DIR / filename
    if not path.exists():
        print(f"MISSING FILE: {path}")
        return None
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)

def get_ids(df, candidates):
    column = find_column(df, candidates)
    if column is None:
        return None, None
    values = df[column].dropna().astype(str).str.strip()
    values = values[values != ""]
    return column, set(values)

def check_table(filename, df):
    print(f"\n--- {filename} ---")
    print(f"Rows: {len(df):,} | Columns: {len(df.columns)}")
    print("Columns:", ", ".join(df.columns))

    duplicate_rows = int(df.duplicated().sum())
    print(f"Exact duplicate rows: {duplicate_rows:,}")

    empty_cells = int(
        (df.isna() | df.apply(
            lambda col: col.astype(str).str.strip().eq("")
        )).sum().sum()
    )
    print(f"Missing/empty cells: {empty_cells:,}")

    id_candidates = [
        "id", "targetId", "target_id", "gene_id", "geneId",
        "diseaseId", "disease_id", "variant_id", "variantId",
        "study_id", "studyId", "rsid"
    ]
    id_col = find_column(df, id_candidates)

    if id_col:
        ids = df[id_col].dropna().astype(str).str.strip()
        ids = ids[ids != ""]
        print(f"Candidate ID column: {id_col}")
        print(f"Unique non-empty IDs: {ids.nunique():,}")
        print(f"Repeated IDs in this column: {ids.duplicated().sum():,}")

def check_relationship(df, source_candidates, valid_ids, label):
    if df is None or valid_ids is None:
        print(f"{label}: SKIPPED (table or reference IDs unavailable)")
        return

    col = find_column(df, source_candidates)
    if col is None:
        print(f"{label}: SKIPPED (could not identify the ID column)")
        return

    values = df[col].dropna().astype(str).str.strip()
    values = values[values != ""]
    unmatched = values[~values.isin(valid_ids)]

    print(
        f"{label}: column '{col}', "
        f"{len(unmatched):,} unmatched values "
        f"out of {len(values):,} non-empty values"
    )
    if len(unmatched):
        print("  Examples:", unmatched.drop_duplicates().head(5).tolist())

def check_raw_gwas_files():
    print("\n" + "=" * 70)
    print("RAW GWAS FILE DUPLICATION CHECK")
    print("=" * 70)

    if not RAW_DIR.exists():
        print("Raw data directory not found:", RAW_DIR)
        return

    files = [
        p for p in RAW_DIR.rglob("*")
        if p.is_file()
        and p.suffix.lower() in {".tsv", ".txt", ".csv"}
        and "gwas" in str(p).lower()
    ]

    # Remove repeated paths, even if folder capitalization differs.
    unique_paths = {}
    for path in files:
        key = str(path.resolve()).casefold()
        unique_paths[key] = path

    hashes = {}
    for path in unique_paths.values():
        digest = hashlib.sha256()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                digest.update(chunk)

        value = digest.hexdigest()
        print(f"{path} | SHA256: {value}")
        hashes.setdefault(value, []).append(str(path))

    duplicates_found = False
    for paths in hashes.values():
        if len(paths) > 1:
            duplicates_found = True
            print("\nIdentical-content files:")
            for path in paths:
                print("  ", path)

    if not duplicates_found:
        print("No separate GWAS files with identical content were found.")

def main():
    print("=" * 70)
    print("CLEANED DATA VALIDATION")
    print("=" * 70)

    tables = {}
    for filename in FILES:
        df = read_tsv(filename)
        if df is not None:
            tables[filename] = df
            check_table(filename, df)

    # Check unique IDs for the main entity tables.
    genes = tables.get("genes.tsv")
    diseases = tables.get("diseases.tsv")
    variants = tables.get("gwas_variants.tsv")
    studies = tables.get("gwas_studies.tsv")

    gene_col, gene_ids = get_ids(genes, ["id", "targetId", "gene_id"]) if genes is not None else (None, None)
    disease_col, disease_ids = get_ids(diseases, ["id", "diseaseId", "disease_id"]) if diseases is not None else (None, None)
    variant_col, variant_ids = get_ids(variants, ["variant_id", "variantId", "id", "variant", "rsid", "SNPS"]) if variants is not None else (None, None)
    study_col, study_ids = get_ids(studies, ["study_id", "studyId", "id", "PUBMEDID"]) if studies is not None else (None, None)

    print("\n" + "=" * 70)
    print("RELATIONSHIP VALIDATION")
    print("=" * 70)

    associations = tables.get("gene_disease_associations.tsv")
    check_relationship(
        associations, ["targetId", "target_id", "gene_id", "geneId"],
        gene_ids, "Association -> Gene"
    )
    check_relationship(
        associations, ["diseaseId", "disease_id"],
        disease_ids, "Association -> Disease"
    )

    hierarchy = tables.get("disease_hierarchy.tsv")
    check_relationship(
        hierarchy, ["diseaseId", "disease_id", "child_id", "childId"],
        disease_ids, "Disease hierarchy -> Child disease"
    )
    check_relationship(
        hierarchy, ["parentId", "parent_id", "parent", "parentDiseaseId"],
        disease_ids, "Disease hierarchy -> Parent disease"
    )

    variant_study = tables.get("gwas_variant_study.tsv")
    check_relationship(
        variant_study, ["variant_id", "variantId", "variant", "rsid", "SNPS"],
        variant_ids, "GWAS variant-study link -> Variant"
    )
    check_relationship(
        variant_study, ["study_id", "studyId", "PUBMEDID", "pubmed_id"],
        study_ids, "GWAS variant-study link -> Study"
    )

    study_disease = tables.get("gwas_study_disease.tsv")
    check_relationship(
        study_disease, ["study_id", "studyId", "PUBMEDID", "pubmed_id"],
        study_ids, "GWAS study-disease link -> Study"
    )

    variant_gene = tables.get("gwas_variant_gene.tsv")
    check_relationship(
        variant_gene, ["variant_id", "variantId", "variant", "rsid", "SNPS"],
        variant_ids, "GWAS variant-gene link -> Variant"
    )

    # GWAS gene links may use gene symbols rather than Open Targets IDs.
    print(
        "Note: GWAS variant-gene links may use gene symbols or other identifiers. "
        "Unmatched gene values need column-level inspection before being treated as errors."
    )

    check_raw_gwas_files()

    print("\n" + "=" * 70)
    print("VALIDATION FINISHED")
    print("Review unmatched IDs and duplicate counts before graph construction.")
    print("=" * 70)

if __name__ == "__main__":
    main()
