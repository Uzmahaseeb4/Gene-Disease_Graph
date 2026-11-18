"""
Gene–Disease Knowledge Graph: data exploration and cleaning
File: 02_explore_clean.py

Expected project layout:
    data/
      raw/
        opentargets/
          target/
          disease/
          association_overall_direct/
          association_by_datatype_direct/   (optional)
        gwas/                               (GWAS Catalog TSV file(s), optional)
      clean/
      reports/

This script reads Open Targets Parquet datasets in batches for the large
association tables to avoid loading the full table into RAM at once.
"""

import os
import re
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
RAW = "data/raw"
CLEAN = "data/clean"
REPORT = "data/reports"

MIN_SCORE = 0.1
GWAS_SIG = 5e-8
BATCH_SIZE = 5000

STD_CHR = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


# ---------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------
def section(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def ensure_dirs():
    for folder in (CLEAN, REPORT):
        os.makedirs(folder, exist_ok=True)


def clean_text(series):
    """Convert a Series to trimmed nullable strings."""
    return series.astype("string").str.strip().replace("", pd.NA)



def profile(df, name):
    """Create a profile report that supports nested Open Targets fields."""
    import json

    os.makedirs(REPORT, exist_ok=True)

    if df is None or df.empty:
        print(f"{name}: empty dataset")
        return

    def safe_nunique(series):
        """Count unique values, including dictionaries and nested lists."""
        def make_hashable(value):
            if isinstance(value, np.ndarray):
                value = value.tolist()

            if isinstance(value, (dict, list, tuple, set)):
                if isinstance(value, set):
                    value = sorted(value, key=str)

                return json.dumps(
                    value,
                    sort_keys=True,
                    default=str
                )

            return value

        try:
            return series.map(make_hashable).nunique(dropna=True)
        except (TypeError, ValueError):
            # Fallback for unusual mixed-type columns.
            return series.map(
                lambda x: str(x) if isinstance(
                    x, (dict, list, tuple, set, np.ndarray)
                ) else x
            ).nunique(dropna=True)

    summary = pd.DataFrame({
        "dtype": df.dtypes.astype(str),
        "null_pct": (df.isna().mean() * 100).round(2),
        "n_unique": pd.Series({
            column: safe_nunique(df[column])
            for column in df.columns
        }),
    })

    path = os.path.join(REPORT, f"profile_{name}.csv")
    summary.to_csv(path, index_label="column")

    print(f"{name}: {len(df):,} rows x {len(df.columns)} cols")
    print(summary.head(12).to_string())
    print(f"Profile saved to: {path}")



def save(df, filename):
    os.makedirs(CLEAN, exist_ok=True)
    path = os.path.join(CLEAN, filename)
    df.to_csv(path, sep="\t", index=False)
    print(f"  -> saved {path}  ({len(df):,} rows)")
    return path


def find_dataset_dir(folder):
    """Find an Open Targets Parquet dataset directory."""
    base = os.path.join(RAW, "opentargets", folder)
    if os.path.isdir(base):
        return base

    # Handle common alternative nesting/layouts.
    candidates = [
        os.path.join(RAW, folder),
        os.path.join(RAW, "opentargets", f"{folder}.parquet"),
    ]
    for candidate in candidates:
        if os.path.isdir(candidate):
            return candidate
    raise FileNotFoundError(
        f"Open Targets dataset folder not found for '{folder}'. "
        f"Expected: {base}"
    )


def read_ot(folder, columns=None):
    """
    Read an Open Targets Parquet dataset.

    Use only for smaller datasets (for example target and disease). Avoid
    using this helper for the full association dataset; use iter_ot_batches.
    """
    path = find_dataset_dir(folder)
    dataset = ds.dataset(path, format="parquet")
    if columns is not None:
        selected = [c for c in columns if c in dataset.schema.names]
    else:
        selected = None
    return dataset.to_table(columns=selected).to_pandas()


def iter_ot_batches(folder, columns=None, batch_size=BATCH_SIZE):
    """Yield smaller DataFrames from an Open Targets Parquet dataset."""
    path = find_dataset_dir(folder)
    dataset = ds.dataset(path, format="parquet")
    if columns is not None:
        selected = [c for c in columns if c in dataset.schema.names]
    else:
        selected = None
    scanner = dataset.scanner(
        columns=selected,
        batch_size=batch_size,
        use_threads=False,
    )
    for batch in scanner.to_batches():
        yield batch.to_pandas()


def first_value(value):
    """Get the first item from a list/array-like value, if possible."""
    if isinstance(value, (list, tuple, np.ndarray)):
        return value[0] if len(value) else None
    return value


def get_first_location(value):
    """Extract the first genomic-location record from nested data."""
    value = first_value(value)
    if isinstance(value, dict):
        return value
    return {}


def parse_chromosome(value):
    if pd.isna(value):
        return pd.NA
    value = str(value).strip()
    value = re.sub(r"^chr", "", value, flags=re.IGNORECASE)
    value = value.upper()
    if value == "M":
        value = "MT"
    return value if value in STD_CHR else pd.NA


def to_numeric(series):
    return pd.to_numeric(series, errors="coerce")


def choose_column(columns, candidates):
    """Return the first candidate present in a list of column names."""
    lookup = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def get_gwas_files():
    """Locate GWAS Catalog TSV/TXT files beneath data/raw."""
    patterns = [
        os.path.join(RAW, "gwas", "**", "*.tsv"),
        os.path.join(RAW, "gwas", "**", "*.tsv.gz"),
        os.path.join(RAW, "gwas", "**", "*.txt"),
        os.path.join(RAW, "GWAS", "**", "*.tsv"),
        os.path.join(RAW, "**", "*gwas*.tsv"),
        os.path.join(RAW, "**", "*GWAS*.tsv"),
    ]
    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern, recursive=True))
    # Exclude generated clean/report files and deduplicate.
    files = sorted({
        os.path.normpath(f) for f in files
        if os.path.isfile(f)
        and os.path.normpath(CLEAN) not in os.path.normpath(f)
        and os.path.normpath(REPORT) not in os.path.normpath(f)
    })
    
    unique_files = []
    seen = set()

    for file_path in files:
        key = os.path.normcase(os.path.abspath(file_path))
        if key not in seen:
            seen.add(key)
            unique_files.append(file_path)

    return unique_files



def read_gwas_file(path):
    """Read a potentially large GWAS TSV with tolerant encoding."""
    compression = "gzip" if path.endswith(".gz") else "infer"
    try:
        return pd.read_csv(
            path, sep="\t", dtype="string", low_memory=False,
            compression=compression,
        )
    except UnicodeDecodeError:
        return pd.read_csv(
            path, sep="\t", dtype="string", low_memory=False,
            encoding="latin-1", compression=compression,
        )


# ---------------------------------------------------------------------
# Section 1: Genes / targets
# ---------------------------------------------------------------------
def process_genes():
    section("1. TARGETS (genes)")
    df = read_ot(
        "target",
        columns=["id", "approvedSymbol", "approvedName", "biotype",
                 "genomicLocation", "genomicLocations"],
    )
    profile(df, "target_raw")

    out = pd.DataFrame()
    out["ensembl_id"] = clean_text(df.get("id", pd.Series(dtype="string")))
    out["gene_symbol"] = clean_text(
        df.get("approvedSymbol", pd.Series(dtype="string"))
    )
    out["gene_name"] = clean_text(
        df.get("approvedName", pd.Series(dtype="string"))
    )
    out["biotype"] = clean_text(
        df.get("biotype", pd.Series(dtype="string"))
    )

    # Open Targets releases may use either genomicLocation or genomicLocations.
    location_col = (
        "genomicLocation" if "genomicLocation" in df.columns else
        "genomicLocations" if "genomicLocations" in df.columns else None
    )
    if location_col:
        locations = df[location_col].apply(get_first_location)
        out["chromosome"] = locations.apply(
            lambda x: parse_chromosome(x.get("chromosome"))
        )
        out["start"] = pd.to_numeric(
            locations.apply(lambda x: x.get("start")), errors="coerce"
        ).astype("Int64")
        out["end"] = pd.to_numeric(
            locations.apply(lambda x: x.get("end")), errors="coerce"
        ).astype("Int64")
        out["strand"] = clean_text(
            locations.apply(lambda x: x.get("strand")).astype("string")
        )
    else:
        out["chromosome"] = pd.NA
        out["start"] = pd.Series(pd.NA, index=out.index, dtype="Int64")
        out["end"] = pd.Series(pd.NA, index=out.index, dtype="Int64")
        out["strand"] = pd.NA

    out = out.dropna(subset=["ensembl_id"]).drop_duplicates("ensembl_id")
    out = out[out["ensembl_id"].str.match(r"^ENSG\d+", na=False)].copy()

    print(f"  kept {len(out):,}/{len(df):,} genes; biotypes:")
    print(out["biotype"].value_counts(dropna=False).head(10).to_string())
    save(out, "genes.tsv")
    return set(out["ensembl_id"].dropna().astype(str))


# ---------------------------------------------------------------------
# Section 2: Diseases and disease edges
# ---------------------------------------------------------------------
def process_diseases():
    section("2. DISEASES")
    df = read_ot(
        "disease",
        columns=["id", "name", "description", "parents", "therapeuticAreas",
                 "ontology"],
    )
    profile(df, "disease_raw")

    diseases = pd.DataFrame({
        "disease_id": clean_text(df.get("id", pd.Series(dtype="string"))),
        "disease_name": clean_text(df.get("name", pd.Series(dtype="string"))),
        "description": clean_text(
            df.get("description", pd.Series(dtype="string"))
        ),
        "ontology": clean_text(df.get("ontology", pd.Series(dtype="string"))),
    })
    diseases = diseases.dropna(subset=["disease_id"])
    diseases = diseases.drop_duplicates("disease_id")
    save(diseases, "diseases.tsv")

    valid_diseases = set(diseases["disease_id"].astype(str))

    hierarchy_rows = []
    area_rows = []
    for _, row in df.iterrows():
        disease_id = row.get("id")
        if pd.isna(disease_id):
            continue

        parents = row.get("parents")
        if isinstance(parents, (list, tuple, np.ndarray)):
            for parent in parents:
                parent_id = first_value(parent)
                if isinstance(parent_id, dict):
                    parent_id = parent_id.get("id")
                if parent_id and str(parent_id) in valid_diseases:
                    hierarchy_rows.append({
                        "child_disease_id": str(disease_id),
                        "parent_disease_id": str(parent_id),
                    })

        areas = row.get("therapeuticAreas")
        if isinstance(areas, (list, tuple, np.ndarray)):
            for area in areas:
                area_id = first_value(area)
                if isinstance(area_id, dict):
                    area_id = area_id.get("id")
                if area_id:
                    area_rows.append({
                        "disease_id": str(disease_id),
                        "therapeutic_area_id": str(area_id),
                    })

    hierarchy = pd.DataFrame(
        hierarchy_rows,
        columns=["child_disease_id", "parent_disease_id"],
    ).drop_duplicates()
    area_edges = pd.DataFrame(
        area_rows,
        columns=["disease_id", "therapeutic_area_id"],
    ).drop_duplicates()

    save(hierarchy, "disease_hierarchy.tsv")
    save(area_edges, "disease_therapeutic_area.tsv")
    return valid_diseases


# ---------------------------------------------------------------------
# Section 3: Gene–disease associations
# ---------------------------------------------------------------------
def process_associations(valid_genes, valid_diseases):
    section("3. GENE-DISEASE ASSOCIATIONS (direct, overall)")

    folder = "association_overall_direct"
    path = find_dataset_dir(folder)
    dataset = ds.dataset(path, format="parquet")
    schema_names = dataset.schema.names
    print(f"Association dataset columns: {schema_names}")

    score_candidates = ["score", "associationScore", "overallScore",
                        "aggregationValue"]
    score_col = choose_column(schema_names, score_candidates)

    # Fallback: use a numeric field that looks like a score.
    if score_col is None:
        for field in dataset.schema:
            if field.name in {"targetId", "diseaseId", "evidenceCount"}:
                continue
            if pa.types.is_floating(field.type):
                score_col = field.name
                break

    if score_col is None:
        raise ValueError(
            "Could not identify the association score column. "
            f"Available columns: {schema_names}"
        )
    print(f"Using '{score_col}' as the association score.")

    required = ["targetId", "diseaseId", score_col]
    missing = [c for c in required if c not in schema_names]
    if missing:
        raise ValueError(f"Association dataset is missing columns: {missing}")

    columns = required + (["evidenceCount"] if "evidenceCount" in schema_names else [])
    scanner = dataset.scanner(
        columns=columns,
        batch_size=BATCH_SIZE,
        use_threads=False,
    )

    association_parts = []
    scanned = 0
    retained = 0
    gene_ids = set(map(str, valid_genes))
    disease_ids = set(map(str, valid_diseases))

    for batch in scanner.to_batches():
        chunk = batch.to_pandas()
        scanned += len(chunk)

        chunk = chunk.rename(columns={
            "targetId": "ensembl_id",
            "diseaseId": "disease_id",
            score_col: "score",
            "evidenceCount": "evidence_count",
        })
        chunk["ensembl_id"] = clean_text(chunk["ensembl_id"])
        chunk["disease_id"] = clean_text(chunk["disease_id"])
        chunk["score"] = to_numeric(chunk["score"])
        chunk = chunk.dropna(subset=["ensembl_id", "disease_id", "score"])
        chunk = chunk[
            chunk["score"].between(0, 1, inclusive="both")
            & (chunk["score"] >= MIN_SCORE)
            & chunk["ensembl_id"].isin(gene_ids)
            & chunk["disease_id"].isin(disease_ids)
        ].copy()

        if not chunk.empty:
            keep_cols = ["ensembl_id", "disease_id", "score"]
            if "evidence_count" in chunk.columns:
                chunk["evidence_count"] = to_numeric(chunk["evidence_count"])
                keep_cols.append("evidence_count")
            association_parts.append(chunk[keep_cols])
            retained += len(chunk)

        if scanned and scanned % 100_000 < BATCH_SIZE:
            print(f"  Processed {scanned:,} association rows...")

    print(f"  Total association rows scanned: {scanned:,}")
    if association_parts:
        associations = pd.concat(association_parts, ignore_index=True)
    else:
        associations = pd.DataFrame(
            columns=["ensembl_id", "disease_id", "score"]
        )

    if not associations.empty:
        # Keep the strongest score for each gene–disease pair.
        agg = {"score": "max"}
        if "evidence_count" in associations.columns:
            agg["evidence_count"] = "max"
        associations = (
            associations.groupby(["ensembl_id", "disease_id"], as_index=False)
            .agg(agg)
        )
        associations["score"] = associations["score"].round(6)

    print(f"  Associations retained after cleaning: {len(associations):,}")
    if not associations.empty:
        print("  Score distribution:")
        print(associations["score"].describe().round(4).to_string())

    save(associations, "gene_disease_associations.tsv")

    # Optional: merge per-datatype association scores.
    datatype_folder = "association_by_datatype_direct"
    try:
        datatype_path = find_dataset_dir(datatype_folder)
        dt_dataset = ds.dataset(datatype_path, format="parquet")
        dt_schema = dt_dataset.schema.names
        print(f"  Merging per-datatype scores from '{datatype_folder}'")

        target_col = choose_column(dt_schema, ["targetId", "ensembl_id"])
        disease_col = choose_column(dt_schema, ["diseaseId", "disease_id"])
        datatype_col = choose_column(
            dt_schema, ["datatypeId", "datatype", "dataTypeId", "datasourceId"]
        )
        value_col = choose_column(
            dt_schema, ["aggregationValue", "score", "associationScore"]
        )

        if not all([target_col, disease_col, datatype_col, value_col]):
            print(
                "  Skipping datatype scores: required target, disease, "
                f"datatype, or score column not found. Columns: {dt_schema}"
            )
        else:
            dt_scanner = dt_dataset.scanner(
                columns=[target_col, disease_col, datatype_col, value_col],
                batch_size=BATCH_SIZE,
                use_threads=False,
            )
            pair_table = associations[["ensembl_id", "disease_id"]].drop_duplicates()
            datatype_parts = []
            dt_scanned = 0

            for batch in dt_scanner.to_batches():
                chunk = batch.to_pandas()
                dt_scanned += len(chunk)
                chunk = chunk.rename(columns={
                    target_col: "ensembl_id",
                    disease_col: "disease_id",
                    datatype_col: "datatype",
                    value_col: "datatype_score",
                })
                chunk["ensembl_id"] = clean_text(chunk["ensembl_id"])
                chunk["disease_id"] = clean_text(chunk["disease_id"])
                chunk["datatype"] = clean_text(chunk["datatype"])
                chunk["datatype_score"] = to_numeric(chunk["datatype_score"])
                chunk = chunk.dropna(
                    subset=["ensembl_id", "disease_id", "datatype",
                            "datatype_score"]
                )
                chunk = chunk[chunk["datatype_score"].between(0, 1)]
                if not chunk.empty:
                    chunk = chunk.merge(
                        pair_table,
                        on=["ensembl_id", "disease_id"],
                        how="inner",
                    )
                    if not chunk.empty:
                        datatype_parts.append(chunk)

                if dt_scanned and dt_scanned % 500_000 < BATCH_SIZE:
                    print(f"    Processed {dt_scanned:,} datatype rows...")

            if datatype_parts:
                dt = pd.concat(datatype_parts, ignore_index=True)
                dt = (
                    dt.groupby(
                        ["ensembl_id", "disease_id", "datatype"],
                        as_index=False,
                    )["datatype_score"].max()
                )
                wide = dt.pivot_table(
                    index=["ensembl_id", "disease_id"],
                    columns="datatype",
                    values="datatype_score",
                    aggfunc="max",
                ).reset_index()
                wide.columns = [
                    str(c) if not isinstance(c, tuple) else "_".join(map(str, c))
                    for c in wide.columns
                ]
                wide = wide.rename(columns={
                    c: f"score_{c}" for c in wide.columns
                    if c not in {"ensembl_id", "disease_id"}
                })
                merged = associations.merge(
                    wide, on=["ensembl_id", "disease_id"], how="left"
                )
                save(merged, "gene_disease_associations.tsv")
                print(f"  Added per-datatype score columns: {len(wide.columns) - 2}")
            else:
                print("  No matching per-datatype score rows found.")
    except FileNotFoundError:
        print("  Optional datatype dataset not found; skipping datatype scores.")


# ---------------------------------------------------------------------
# Section 4: GWAS Catalog
# ---------------------------------------------------------------------
def process_gwas(valid_genes, valid_diseases):
    section("4. GWAS CATALOG")

    files = get_gwas_files()
    if not files:
        print(
            "No GWAS Catalog TSV/TXT file found under data/raw/gwas or "
            "data/raw. Skipping GWAS processing."
        )
        return

    print("GWAS file(s) found:")
    for f in files:
        print(f"  - {f}")

    frames = []
    for path in files:
        try:
            frame = read_gwas_file(path)
            frame["__source_file"] = os.path.basename(path)
            frames.append(frame)
            print(f"  Loaded {os.path.basename(path)}: {len(frame):,} rows")
        except Exception as exc:
            warnings.warn(f"Could not read GWAS file {path}: {exc}")

    if not frames:
        print("No GWAS files could be read; skipping GWAS processing.")
        return

    gwas = pd.concat(frames, ignore_index=True, sort=False)
    profile(gwas, "gwas_raw")

    # Normalise column lookup while preserving original names.
    cols = list(gwas.columns)
    variant_col = choose_column(cols, [
        "SNPS", "SNP_ID_CURRENT", "SNP_ID", "VARIANT_ID", "CHR_POS",
        "CHR_ID", "variant_id", "variant",
    ])
    study_col = choose_column(cols, [
        "STUDY ACCESSION", "STUDY_ACCESSION", "GCST", "study_accession",
        "STUDY",
    ])
    trait_col = choose_column(cols, [
        "DISEASE/TRAIT", "DISEASE_TRAIT", "MAPPED_TRAIT", "TRAIT",
        "disease_trait",
    ])
    pvalue_col = choose_column(cols, ["P-VALUE", "P_VALUE", "PVALUE", "p_value"])
    author_col = choose_column(cols, ["FIRST AUTHOR", "FIRST_AUTHOR", "AUTHOR"])
    date_col = choose_column(cols, ["DATE", "PUBLICATION DATE", "PUBMED DATE"])
    gene_col = choose_column(cols, [
        "MAPPED_GENE", "MAPPED GENE", "REPORTED GENE(S)", "REPORTED_GENE",
        "GENE", "mapped_gene",
    ])
    chr_col = choose_column(cols, ["CHR_ID", "CHR", "CHROMOSOME"])
    pos_col = choose_column(cols, ["CHR_POS", "POSITION", "POS", "CHROMOSOMAL_POSITION"])
    mapped_uri_col = choose_column(cols, ["MAPPED_TRAIT_URI", "MAPPED_TRAIT_URI"])
    pubmed_col = choose_column(cols, ["PUBMEDID", "PUBMED_ID", "PUBMED ID"])

    if study_col is None:
        # If the file has no study accession, use the source filename plus row ID.
        gwas["__study_id"] = (
            gwas["__source_file"].astype(str) + "_row_" +
            pd.Series(np.arange(len(gwas)), index=gwas.index).astype(str)
        )
        study_col = "__study_id"

    if variant_col is None:
        if chr_col and pos_col:
            gwas["__variant_id"] = (
                gwas[chr_col].astype("string").str.strip() + ":" +
                gwas[pos_col].astype("string").str.strip()
            )
            variant_col = "__variant_id"
        else:
            gwas["__variant_id"] = pd.NA
            variant_col = "__variant_id"

    # Variant table.
    variants = pd.DataFrame({
        "variant_id": clean_text(gwas[variant_col]),
    })
    if chr_col:
        variants["chromosome"] = gwas[chr_col].apply(parse_chromosome)
    else:
        variants["chromosome"] = pd.NA
    if pos_col:
        variants["position"] = pd.to_numeric(gwas[pos_col], errors="coerce").astype("Int64")
    else:
        variants["position"] = pd.Series(pd.NA, index=gwas.index, dtype="Int64")
    variants = variants.dropna(subset=["variant_id"]).drop_duplicates("variant_id")
    save(variants, "gwas_variants.tsv")

    # Study table.
    studies = pd.DataFrame({"study_id": clean_text(gwas[study_col])})
    studies["first_author"] = clean_text(gwas[author_col]) if author_col else pd.NA
    studies["publication_date"] = clean_text(gwas[date_col]) if date_col else pd.NA
    studies["pubmed_id"] = clean_text(gwas[pubmed_col]) if pubmed_col else pd.NA
    studies = studies.dropna(subset=["study_id"]).drop_duplicates("study_id")
    save(studies, "gwas_studies.tsv")

    # Variant–study edges and GWAS study–trait/disease edges.
    pvalues = (
        pd.to_numeric(
            gwas[pvalue_col].astype("string").str.replace(",", "", regex=False),
            errors="coerce",
        ) if pvalue_col else pd.Series(np.nan, index=gwas.index)
    )

    variant_study = pd.DataFrame({
        "variant_id": clean_text(gwas[variant_col]),
        "study_id": clean_text(gwas[study_col]),
        "p_value": pvalues,
    })
    if pvalue_col:
        variant_study["is_significant"] = variant_study["p_value"] <= GWAS_SIG
    else:
        variant_study["is_significant"] = pd.NA
    variant_study = variant_study.dropna(
        subset=["variant_id", "study_id"]
    ).drop_duplicates(["variant_id", "study_id"])
    save(variant_study, "gwas_variant_study.tsv")

    study_disease = pd.DataFrame({
        "study_id": clean_text(gwas[study_col]),
        "trait": clean_text(gwas[trait_col]) if trait_col else pd.NA,
        "mapped_trait_uri": clean_text(gwas[mapped_uri_col]) if mapped_uri_col else pd.NA,
    })
    study_disease = study_disease.dropna(subset=["study_id"])
    study_disease = study_disease.drop_duplicates(
        ["study_id", "trait", "mapped_trait_uri"]
    )
    save(study_disease, "gwas_study_disease.tsv")

    # Map GWAS gene symbols to Open Targets gene IDs where possible.
    # The raw GWAS Catalog often supplies gene symbols, not Ensembl IDs.
    gene_rows = []
    if gene_col:
        gene_values = clean_text(gwas[gene_col])
        for idx, value in gene_values.dropna().items():
            # GWAS values can contain multiple genes separated by commas,
            # semicolons, slashes, or " - ".
            tokens = re.split(r"[;,/|]", str(value))
            variant_id = gwas.at[idx, variant_col]
            study_id = gwas.at[idx, study_col]
            for token in tokens:
                token = token.strip()
                if token:
                    gene_rows.append({
                        "variant_id": variant_id,
                        "study_id": study_id,
                        "gene_symbol": token,
                    })

    variant_gene = pd.DataFrame(
        gene_rows, columns=["variant_id", "study_id", "gene_symbol"]
    ).dropna(subset=["variant_id", "gene_symbol"]).drop_duplicates()
    save(variant_gene, "gwas_variant_gene.tsv")

    print("\nGWAS summary:")
    print(f"  Unique variants: {variants['variant_id'].nunique():,}")
    print(f"  Unique studies: {studies['study_id'].nunique():,}")
    print(f"  Variant-study links: {len(variant_study):,}")
    print(f"  Study-trait links: {len(study_disease):,}")
    print(f"  Variant-gene links: {len(variant_gene):,}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ensure_dirs()

    genes = process_genes()
    diseases = process_diseases()
    process_associations(genes, diseases)
    process_gwas(genes, diseases)

    section("CLEANING COMPLETE")
    print(f"Clean files:   {os.path.abspath(CLEAN)}")
    print(f"Report files:  {os.path.abspath(REPORT)}")
    print("Review the generated TSV files and profile reports before graph import.")


if __name__ == "__main__":
    main()
