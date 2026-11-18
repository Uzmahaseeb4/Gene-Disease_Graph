"""
Build a Gene-Disease Knowledge Graph in Neo4j.

Input:
    data/clean/genes.tsv
    data/clean/diseases.tsv
    data/clean/gene_disease_associations.tsv
    data/clean/gwas_variants.tsv
    data/clean/gwas_studies.tsv
    data/clean/gwas_variant_study.tsv
    data/clean/gwas_study_disease.tsv
    data/clean/gwas_variant_gene.tsv
    data/clean/disease_hierarchy.tsv

Run:
    python neoj4_fixed.py
"""

import os
import re
from pathlib import Path

import pandas as pd
from neo4j import GraphDatabase

# ---------------------------------------------------------
# 1. CONFIGURATION — update these for your Neo4j database
# ---------------------------------------------------------

NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "project_big_data")

# Set to the database name you created.
NEO4J_DATABASE = os.getenv("NEO4J_DATABASE", "neo4j")

CLEAN_DIR = Path("data/clean")
BATCH_SIZE = 2000


# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------

def read_table(filename):
    path = CLEAN_DIR / filename

    if not path.exists():
        print(f"WARNING: Missing file: {path}")
        return None

    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df.columns = [str(c).strip() for c in df.columns]

    print(f"Loaded {filename}: {len(df):,} rows")
    print("  Columns:", df.columns.tolist())
    return df


def read_first_available(filenames):
    """Read the first existing file from a list of candidate filenames."""
    for filename in filenames:
        if (CLEAN_DIR / filename).exists():
            if filename != filenames[0]:
                print(
                    f"NOTE: {filenames[0]} not found; "
                    f"using {filename} instead."
                )
            return read_table(filename)

    print(f"WARNING: None of these files were found: {filenames}")
    return None


def normalize_name(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def find_column(df, candidates, required=True):
    """Find a column using case- and punctuation-insensitive aliases."""
    if df is None:
        return None

    lookup = {normalize_name(col): col for col in df.columns}

    for candidate in candidates:
        result = lookup.get(normalize_name(candidate))
        if result is not None:
            return result

    if required:
        raise ValueError(
            f"Could not find any of {candidates}. "
            f"Available columns: {df.columns.tolist()}"
        )

    return None


def clean_value(value):
    """Convert a DataFrame value to a Neo4j-compatible property value."""
    if value is None or pd.isna(value):
        return None

    value = str(value).strip()

    if value == "" or value.lower() in {"nan", "none", "<na>"}:
        return None

    return value


def make_properties(row):
    """Keep useful source columns as string properties."""
    properties = {}

    for col, value in row.items():
        value = clean_value(value)

        if value is not None:
            properties[col] = value

    return properties


def batches(items, size=BATCH_SIZE):
    for start in range(0, len(items), size):
        yield items[start:start + size]


def rows_with_keys(df, key_column, extra_properties=None):
    """Convert rows into records with a required graph identifier."""
    if df is None:
        return []

    extra_properties = extra_properties or {}
    records = []

    for _, row in df.iterrows():
        key = clean_value(row.get(key_column))

        if key is None:
            continue

        props = make_properties(row)

        for new_name, source_column in extra_properties.items():
            value = clean_value(row.get(source_column))
            if value is not None:
                props[new_name] = value

        records.append({"key": key, "props": props})

    return records


def create_nodes(session, label, key_property, records):
    """Create or update nodes in batches. Labels are fixed in this script."""
    if not records:
        print(f"{label}: no valid node records")
        return

    query = f"""
    UNWIND $rows AS row
    MERGE (n:{label} {{{key_property}: row.key}})
    SET n += row.props
    """

    total = 0
    for batch in batches(records):
        session.run(query, rows=batch).consume()
        total += len(batch)

    print(f"{label} nodes processed: {total:,}")


def relationship_rows(df, left_candidates, right_candidates,
                      left_name, right_name, extra_candidates=None):
    """Prepare source/target IDs and optional relationship properties."""
    if df is None:
        return []

    left_col = find_column(df, left_candidates)
    right_col = find_column(df, right_candidates)

    extra_col = None
    if extra_candidates:
        extra_col = find_column(df, extra_candidates, required=False)

    records = []

    for _, row in df.iterrows():
        left_id = clean_value(row.get(left_col))
        right_id = clean_value(row.get(right_col))

        if left_id is None or right_id is None:
            continue

        item = {left_name: left_id, right_name: right_id}

        if extra_col:
            extra_value = clean_value(row.get(extra_col))
            if extra_value is not None:
                item["extra"] = extra_value

        records.append(item)

    return records


def create_relationships(session, description, query, records):
    if not records:
        print(f"{description}: no valid relationship rows")
        return

    total = 0

    for batch in batches(records):
        result = session.run(query, rows=batch)
        summary = result.consume()
        total += summary.counters.relationships_created
        total += summary.counters.properties_set * 0  # Keep count creation-based.

    print(
        f"{description}: {len(records):,} rows submitted; "
        f"{total:,} new relationships created"
    )


# ---------------------------------------------------------
# 3. LOAD CLEANED TABLES
# ---------------------------------------------------------

def main():
    genes = read_table("genes.tsv")
    diseases = read_first_available(["diseases.tsv", "disease_therapeutic_area.tsv"])
    associations = read_table("gene_disease_associations.tsv")

    variants = read_table("gwas_variants.tsv")
    studies = read_table("gwas_studies.tsv")
    variant_study = read_table("gwas_variant_study.tsv")
    study_disease = read_table("gwas_study_disease.tsv")
    variant_gene = read_table("gwas_variant_gene.tsv")
    hierarchy = read_table("disease_hierarchy.tsv")

    if genes is None or diseases is None:
        raise FileNotFoundError(
            "genes.tsv and a disease table (diseases.tsv or "
            "disease_therapeutic_area.tsv) are required. "
            "Check data/clean/ and run the cleaning script if needed."
        )

    # Identify the primary IDs from the actual cleaned column names.
    gene_id_col = find_column(
        genes,
        ["ensembl_id", "id", "targetId", "gene_id", "geneId"]
    )
    disease_id_col = find_column(
        diseases, ["id", "diseaseId", "disease_id"]
    )

    gene_symbol_col = find_column(
        genes, ["approvedSymbol", "symbol", "gene_symbol"],
        required=False
    )

    # Build Gene nodes. Also store a standard 'symbol' property for matching.
    gene_extra = {}
    if gene_symbol_col:
        gene_extra["symbol"] = gene_symbol_col

    gene_records = rows_with_keys(
        genes, gene_id_col, extra_properties=gene_extra
    )

    disease_records = rows_with_keys(diseases, disease_id_col)

    # Identify variant and study identifiers.
    variant_id_col = (
        find_column(
            variants,
            ["variant_id", "variantId", "id", "variant", "rsid", "SNPS"],
            required=False
        ) if variants is not None else None
    )

    study_id_col = (
        find_column(
            studies,
            ["study_id", "studyId", "id", "study", "PUBMEDID"],
            required=False
        ) if studies is not None else None
    )

    variant_records = (
        rows_with_keys(variants, variant_id_col)
        if variants is not None and variant_id_col else []
    )

    study_records = (
        rows_with_keys(studies, study_id_col)
        if studies is not None and study_id_col else []
    )

    if variants is not None and not variant_id_col:
        print("WARNING: Could not identify the primary ID in gwas_variants.tsv.")
    if studies is not None and not study_id_col:
        print("WARNING: Could not identify the primary ID in gwas_studies.tsv.")

    # -----------------------------------------------------
    # 4. CONNECT TO NEO4J AND CREATE CONSTRAINTS & INDEXES
    # -----------------------------------------------------

    if NEO4J_PASSWORD == "CHANGE_ME":
        raise ValueError(
            "Set your Neo4j password in NEO4J_PASSWORD near the top "
            "of this script before running it."
        )

    driver = GraphDatabase.driver(
        NEO4J_URI,
        auth=(NEO4J_USER, NEO4J_PASSWORD)
    )

    try:
        driver.verify_connectivity()
        print("\nConnected to Neo4j successfully.")

        with driver.session(database=NEO4J_DATABASE) as session:

            constraints_and_indexes = [
                "CREATE CONSTRAINT gene_id_unique IF NOT EXISTS "
                "FOR (n:Gene) REQUIRE n.gene_id IS UNIQUE",

                "CREATE CONSTRAINT disease_id_unique IF NOT EXISTS "
                "FOR (n:Disease) REQUIRE n.disease_id IS UNIQUE",

                "CREATE CONSTRAINT variant_id_unique IF NOT EXISTS "
                "FOR (n:Variant) REQUIRE n.variant_id IS UNIQUE",

                "CREATE CONSTRAINT study_id_unique IF NOT EXISTS "
                "FOR (n:Study) REQUIRE n.study_id IS UNIQUE",

                "CREATE INDEX gene_symbol_idx IF NOT EXISTS "
                "FOR (n:Gene) ON (n.symbol)",
            ]

            for query in constraints_and_indexes:
                session.run(query).consume()

            print("\nCreating nodes...")
            create_nodes(
                session, "Gene", "gene_id", gene_records
            )
            create_nodes(
                session, "Disease", "disease_id", disease_records
            )
            create_nodes(
                session, "Variant", "variant_id", variant_records
            )
            create_nodes(
                session, "Study", "study_id", study_records
            )

            # -------------------------------------------------
            # 5. GENE-DISEASE ASSOCIATIONS
            # -------------------------------------------------

            assoc_records = relationship_rows(
                associations,
                # This cleaned file uses "ensembl_id" for the gene key.
                ["ensembl_id", "targetId", "target_id", "gene_id", "geneId"],
                ["diseaseId", "disease_id"],
                "gene_id", "disease_id",
                ["associationScore", "score", "association_score"]
            )

            assoc_query = """
            UNWIND $rows AS row
            MATCH (g:Gene {gene_id: row.gene_id})
            MATCH (d:Disease {disease_id: row.disease_id})
            MERGE (g)-[r:HAS_ASSOCIATION]->(d)
            FOREACH (_ IN CASE WHEN row.extra IS NULL THEN [] ELSE [1] END |
                SET r.score = toFloat(row.extra)
            )
            """

            create_relationships(
                session, "Gene-Disease associations",
                assoc_query, assoc_records
            )

            # -------------------------------------------------
            # 6. VARIANT-STUDY RELATIONSHIPS
            # -------------------------------------------------

            if variant_study is not None:
                vs_records = relationship_rows(
                    variant_study,
                    ["variant_id", "variantId", "variant", "rsid", "SNPS"],
                    ["study_id", "studyId", "study", "PUBMEDID", "pubmed_id"],
                    "variant_id", "study_id"
                )

                vs_query = """
                UNWIND $rows AS row
                MATCH (v:Variant {variant_id: row.variant_id})
                MATCH (s:Study {study_id: row.study_id})
                MERGE (v)-[:REPORTED_IN]->(s)
                """

                create_relationships(
                    session, "Variant-Study links", vs_query, vs_records
                )

            # -------------------------------------------------
            # 7. STUDY-DISEASE RELATIONSHIPS (INDEX-OPTIMIZED)
            # -------------------------------------------------

            if study_disease is not None:
                study_col = find_column(
                    study_disease,
                    ["study_id", "studyId", "study", "PUBMEDID", "pubmed_id"],
                )
                uri_col = find_column(
                    study_disease,
                    ["disease_id", "diseaseId", "trait_id", "mapped_trait_uri"],
                    required=False,
                )
                trait_col = find_column(
                    study_disease,
                    ["trait", "disease", "disease_name", "name"],
                    required=False,
                )

                sd_id_records = []
                sd_name_records = []

                for _, row in study_disease.iterrows():
                    s_id = clean_value(row.get(study_col))
                    if not s_id:
                        continue

                    raw_uri = clean_value(row.get(uri_col)) if uri_col else None
                    trait_name = clean_value(row.get(trait_col)) if trait_col else None

                    if raw_uri:
                        for part in raw_uri.split(","):
                            part = part.strip()
                            if part:
                                short_id = part.rstrip("/").split("/")[-1]
                                sd_id_records.append({
                                    "study_id": s_id,
                                    "disease_id": short_id,
                                })
                    elif trait_name:
                        sd_name_records.append({
                            "study_id": s_id,
                            "trait": trait_name.lower(),
                        })

                sd_id_query = """
                UNWIND $rows AS row
                MATCH (s:Study {study_id: row.study_id})
                MATCH (d:Disease {disease_id: row.disease_id})
                MERGE (s)-[:INVESTIGATES]->(d)
                """

                create_relationships(
                    session, "Study-Disease links (by ID)", sd_id_query, sd_id_records
                )

                if sd_name_records:
                    sd_name_query = """
                    UNWIND $rows AS row
                    MATCH (s:Study {study_id: row.study_id})
                    MATCH (d:Disease)
                    WHERE toLower(coalesce(d.name, d.disease_name, d.trait, "")) = row.trait
                    MERGE (s)-[:INVESTIGATES]->(d)
                    """

                    create_relationships(
                        session, "Study-Disease links (by name)", sd_name_query, sd_name_records
                    )

            # -------------------------------------------------
            # 8. VARIANT-GENE RELATIONSHIPS (INDEX-OPTIMIZED)
            # -------------------------------------------------

            if variant_gene is not None:
                vg_records = relationship_rows(
                    variant_gene,
                    ["variant_id", "variantId", "variant", "rsid", "SNPS"],
                    [
                        "gene_id", "geneId", "targetId", "target_id",
                        "gene_symbol", "symbol", "gene", "mapped_gene",
                        "MAPPED_GENE"
                    ],
                    "variant_id", "gene_key"
                )

                vg_id_records = [r for r in vg_records if r["gene_key"].startswith("ENSG")]
                vg_sym_records = [r for r in vg_records if not r["gene_key"].startswith("ENSG")]

                vg_id_query = """
                UNWIND $rows AS row
                MATCH (v:Variant {variant_id: row.variant_id})
                MATCH (g:Gene {gene_id: row.gene_key})
                MERGE (v)-[:ASSOCIATED_WITH_GENE]->(g)
                """

                vg_sym_query = """
                UNWIND $rows AS row
                MATCH (v:Variant {variant_id: row.variant_id})
                MATCH (g:Gene {symbol: row.gene_key})
                MERGE (v)-[:ASSOCIATED_WITH_GENE]->(g)
                """

                if vg_id_records:
                    create_relationships(
                        session, "Variant-Gene links (by ID)", vg_id_query, vg_id_records
                    )
                if vg_sym_records:
                    create_relationships(
                        session, "Variant-Gene links (by symbol)", vg_sym_query, vg_sym_records
                    )

            # -------------------------------------------------
            # 9. DISEASE HIERARCHY
            # -------------------------------------------------

            if hierarchy is not None:
                # The cleaned hierarchy file uses child_disease_id and
                # parent_disease_id. Resolve aliases without assuming a schema.
                parent_col = find_column(
                    hierarchy,
                    [
                        "parent_disease_id",
                        "parentId",
                        "parent_id",
                        "parent",
                        "parentDiseaseId",
                    ],
                    required=False,
                )

                child_col = find_column(
                    hierarchy,
                    [
                        "child_disease_id",
                        "childId",
                        "child_id",
                        "diseaseId",
                        "disease_id",
                    ],
                    required=False,
                )

                if parent_col and child_col:
                    dh_records = relationship_rows(
                        hierarchy,
                        [child_col],
                        [parent_col],
                        "child_id",
                        "parent_id",
                    )

                    dh_query = """
                    UNWIND $rows AS row
                    MATCH (child:Disease {disease_id: row.child_id})
                    MATCH (parent:Disease {disease_id: row.parent_id})
                    MERGE (child)-[:IS_A]->(parent)
                    """

                    create_relationships(
                        session,
                        "Disease hierarchy links",
                        dh_query,
                        dh_records,
                    )
                else:
                    print(
                        "Disease hierarchy skipped: could not identify "
                        "parent and child ID columns."
                    )

            # -------------------------------------------------
            # 10. GRAPH SUMMARY
            # -------------------------------------------------

            print("\nGraph summary:")
            summary_query = """
            MATCH (n)
            RETURN labels(n) AS node_labels, count(*) AS total
            ORDER BY node_labels
            """

            for record in session.run(summary_query):
                print(f"  Nodes {record['node_labels']}: {record['total']:,}")

            rel_query = """
            MATCH ()-[r]->()
            RETURN type(r) AS relationship_type, count(*) AS total
            ORDER BY relationship_type
            """

            for record in session.run(rel_query):
                print(
                    f"  Relationships {record['relationship_type']}: "
                    f"{record['total']:,}"
                )

        print("\nKnowledge graph build finished.")

    finally:
        driver.close()


if __name__ == "__main__":
    main()