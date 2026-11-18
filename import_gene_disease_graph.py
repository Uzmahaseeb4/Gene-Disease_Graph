
import csv
import os
import math
from pathlib import Path

from neo4j import GraphDatabase

# --------------------------------------------------
# CONFIGURATION
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent
IMPORT_DIR = ROOT / "data" / "graph" / "neo4j_import"

URI = os.getenv("NEO4J_URI", "neo4j://127.0.0.1:7687")
USERNAME = os.getenv("NEO4J_USERNAME", "neo4j")
PASSWORD = os.getenv("NEO4J_PASSWORD")

BATCH_SIZE = 5000


def clean_value(value):
    """Convert CSV values to Neo4j-compatible Python values."""
    if value is None or value.strip() == "":
        return None

    value = value.strip()

    try:
        number = float(value)
        if not math.isfinite(number):
            return None
    except ValueError:
        return value

    return value


def parse_number(value, kind):
    value = clean_value(value)

    if value is None:
        return None

    try:
        if kind == "long" or kind == "int":
            return int(value)
        if kind == "double":
            return float(value)
    except (ValueError, OverflowError):
        return None

    return value


def batches(reader, size=BATCH_SIZE):
    batch = []

    for row in reader:
        batch.append(row)

        if len(batch) >= size:
            yield batch
            batch = []

    if batch:
        yield batch


def read_csv_batches(filename):
    path = IMPORT_DIR / filename

    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    handle = path.open(
        "r", encoding="utf-8-sig", newline=""
    )
    reader = csv.DictReader(handle)

    return handle, reader


# --------------------------------------------------
# NODE IMPORTS
# --------------------------------------------------

def import_genes(session):
    handle, reader = read_csv_batches("genes.csv")
    total = 0

    query = """
    UNWIND $rows AS row
    MERGE (g:Gene {gene_id: row.gene_id})
    SET g.symbol = row.symbol,
        g.name = row.name,
        g.biotype = row.biotype,
        g.chromosome = row.chromosome,
        g.start = row.start,
        g.end = row.end,
        g.strand = row.strand
    """

    try:
        for batch in batches(reader):
            rows = []

            for r in batch:
                gene_id = clean_value(
                    r.get("gene_id:ID(Gene)")
                )

                if not gene_id:
                    continue

                rows.append({
                    "gene_id": gene_id,
                    "symbol": clean_value(r.get("symbol")),
                    "name": clean_value(r.get("name")),
                    "biotype": clean_value(r.get("biotype")),
                    "chromosome": clean_value(
                        r.get("chromosome")
                    ),
                    "start": parse_number(
                        r.get("start:long"), "long"
                    ),
                    "end": parse_number(
                        r.get("end:long"), "long"
                    ),
                    "strand": parse_number(
                        r.get("strand"), "int"
                    ),
                })

            if rows:
                session.run(query, rows=rows).consume()
                total += len(rows)
                print(f"Genes processed: {total:,}")

    finally:
        handle.close()

    print(f"Gene import complete: {total:,}")


def import_diseases(session):
    handle, reader = read_csv_batches("diseases.csv")
    total = 0

    query = """
    UNWIND $rows AS row
    MERGE (d:Disease {disease_id: row.disease_id})
    SET d.name = row.name,
        d.description = row.description,
        d.therapeutic_areas = row.therapeutic_areas,
        d.ontology = row.ontology
    """

    try:
        for batch in batches(reader):
            rows = []

            for r in batch:
                disease_id = clean_value(
                    r.get("disease_id:ID(Disease)")
                )

                if not disease_id:
                    continue

                rows.append({
                    "disease_id": disease_id,
                    "name": clean_value(r.get("name")),
                    "description": clean_value(
                        r.get("description")
                    ),
                    "therapeutic_areas": clean_value(
                        r.get("therapeutic_areas")
                    ),
                    "ontology": clean_value(r.get("ontology")),
                })

            if rows:
                session.run(query, rows=rows).consume()
                total += len(rows)
                print(f"Diseases processed: {total:,}")

    finally:
        handle.close()

    print(f"Disease import complete: {total:,}")


# --------------------------------------------------
# ASSOCIATION IMPORT
# --------------------------------------------------

def import_associations(session):
    handle, reader = read_csv_batches("associations.csv")
    total = 0
    skipped = 0

    headers = reader.fieldnames or []

    # Identify typed property columns from the CSV header.
    property_columns = [
        h for h in headers
        if h not in {
            ":START_ID(Gene)",
            ":END_ID(Disease)",
            ":TYPE",
        }
    ]

    query = """
    UNWIND $rows AS row
    MATCH (g:Gene {gene_id: row.gene_id})
    MATCH (d:Disease {disease_id: row.disease_id})
    MERGE (g)-[r:ASSOCIATED_WITH]->(d)
    SET r += row.props
    """

    try:
        for batch in batches(reader):
            rows = []

            for raw in batch:
                gene_id = clean_value(
                    raw.get(":START_ID(Gene)")
                )
                disease_id = clean_value(
                    raw.get(":END_ID(Disease)")
                )

                if not gene_id or not disease_id:
                    skipped += 1
                    continue

                props = {}

                for column in property_columns:
                    # Remove the Neo4j type annotation.
                    key = column.split(":", 1)[0]
                    value = raw.get(column)

                    if column.endswith(":double"):
                        value = parse_number(value, "double")
                    elif column.endswith(":long"):
                        value = parse_number(value, "long")
                    elif column.endswith(":int"):
                        value = parse_number(value, "int")
                    else:
                        # Keep source metadata as a string.
                        # Do not guess how array-like values are encoded.
                        value = clean_value(value)

                    props[key] = value

                rows.append({
                    "gene_id": gene_id,
                    "disease_id": disease_id,
                    "props": props,
                })

            if rows:
                result = session.run(query, rows=rows)
                result.consume()
                total += len(rows)
                print(
                    f"Association rows processed: {total:,}"
                )

    finally:
        handle.close()

    print(f"Association import complete: {total:,}")
    print(f"Rows skipped for missing IDs: {skipped:,}")


# --------------------------------------------------
# VERIFICATION
# --------------------------------------------------

def verify_graph(session):
    query = """
    MATCH (g:Gene)
    WITH count(g) AS genes
    MATCH (d:Disease)
    WITH genes, count(d) AS diseases
    MATCH ()-[r:ASSOCIATED_WITH]->()
    RETURN genes, diseases, count(r) AS associations
    """

    record = session.run(query).single()

    print("\n========== GRAPH SUMMARY ==========")
    print(f"Genes:        {record['genes']:,}")
    print(f"Diseases:     {record['diseases']:,}")
    print(f"Associations: {record['associations']:,}")


def main():
    if not PASSWORD:
        raise RuntimeError(
            "Set NEO4J_PASSWORD in your terminal first."
        )

    required_files = [
        "genes.csv",
        "diseases.csv",
        "associations.csv",
    ]

    for filename in required_files:
        if not (IMPORT_DIR / filename).exists():
            raise FileNotFoundError(
                f"Required import file not found: {filename}"
            )

    print("Connecting to Neo4j...")

    with GraphDatabase.driver(
        URI, auth=(USERNAME, PASSWORD)
    ) as driver:
        driver.verify_connectivity()
        print("Neo4j connection successful.")

        # Create uniqueness constraints before loading nodes.
        with driver.session() as session:
            session.run("""
                CREATE CONSTRAINT gene_id_unique IF NOT EXISTS
                FOR (g:Gene) REQUIRE g.gene_id IS UNIQUE
            """).consume()

            session.run("""
                CREATE CONSTRAINT disease_id_unique IF NOT EXISTS
                FOR (d:Disease) REQUIRE d.disease_id IS UNIQUE
            """).consume()

            print("\nImporting genes...")
            import_genes(session)

            print("\nImporting diseases...")
            import_diseases(session)

            print("\nImporting associations...")
            import_associations(session)

            verify_graph(session)

    print("\nImport pipeline finished.")


if __name__ == "__main__":
    main()
