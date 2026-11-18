import os
from neo4j.exceptions import Neo4jError
from contextlib import asynccontextmanager
from pathlib import Path

from dotenv import load_dotenv
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from neo4j import GraphDatabase


# --------------------------------------------------
# Configuration
# --------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
load_dotenv(BASE_DIR / ".env")

NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD")
NEO4J_DATABASE = os.getenv("NEO4J_DATABASE", "neo4j")

driver = None


# --------------------------------------------------
# Application startup and shutdown
# --------------------------------------------------

@asynccontextmanager
async def lifespan(app: FastAPI):
    global driver

    if not NEO4J_PASSWORD:
        raise RuntimeError(
            "NEO4J_PASSWORD is missing from app/.env"
        )

    driver = GraphDatabase.driver(
        NEO4J_URI,
        auth=(NEO4J_USER, NEO4J_PASSWORD)
    )

    try:
        driver.verify_connectivity()
        print("Connected to Neo4j successfully!")
        yield
    finally:
        driver.close()
        driver = None


app = FastAPI(
    title="Gene-Disease Knowledge Graph API",
    description=(
        "Search genes, diseases, and their associations "
        "stored in Neo4j."
    ),
    version="1.0.0",
    lifespan=lifespan
)

app.mount(
    "/static",
    StaticFiles(directory=BASE_DIR / "static"),
    name="static"
)


# --------------------------------------------------
# Helper function
# --------------------------------------------------

def run_query(query, parameters=None):
    if driver is None:
        raise HTTPException(
            status_code=503,
            detail="Neo4j database is not connected."
        )

    try:
        with driver.session(database=NEO4J_DATABASE) as session:
            result = session.run(query, parameters or {})
            return [record.data() for record in result]

    except Exception as exc:
        print("Neo4j query error:", str(exc))

        raise HTTPException(
            status_code=500,
            detail=(
                "Database query failed. "
                "Check the API terminal for details."
            )
        ) from exc


# --------------------------------------------------
# Home page
# --------------------------------------------------

@app.get("/", include_in_schema=False)
def home():
    return FileResponse(BASE_DIR / "static" / "index.html")


# --------------------------------------------------
# 1. Health endpoint
# --------------------------------------------------

@app.get("/api/health")
def health():
    if driver is None:
        raise HTTPException(
            status_code=503,
            detail="Neo4j is not connected."
        )

    try:
        driver.verify_connectivity()

        return {
            "status": "connected",
            "database": NEO4J_DATABASE,
            "message": "API and Neo4j are working."
        }

    except Exception as exc:
        raise HTTPException(
            status_code=503,
            detail="Could not connect to Neo4j."
        ) from exc


# --------------------------------------------------
# 2. Graph statistics endpoint
# --------------------------------------------------

@app.get("/api/stats")
def graph_stats():

    node_query = """
    MATCH (n)
    UNWIND labels(n) AS label
    RETURN label, count(*) AS total
    ORDER BY total DESC
    """

    relationship_query = """
    MATCH ()-[r]->()
    RETURN type(r) AS relationship_type,
           count(*) AS total
    ORDER BY total DESC
    """

    nodes = run_query(node_query)
    relationships = run_query(relationship_query)

    return {
        "nodes": nodes,
        "relationships": relationships,
        "total_nodes": sum(item["total"] for item in nodes),
        "total_relationships": sum(
            item["total"] for item in relationships
        )
    }


# --------------------------------------------------
# 3. Search endpoint
# --------------------------------------------------

@app.get("/api/search")
def search_graph(
    q: str = Query(..., min_length=1, max_length=100),
    limit: int = Query(20, ge=1, le=100)
):
    search_term = q.strip().lower()

    if not search_term:
        raise HTTPException(
            status_code=422,
            detail="Enter a gene symbol, gene ID, or disease name."
        )

    parameters = {
        "q": search_term,
        "limit": limit
    }

    # Search genes by symbol, name, or identifier.
    gene_query = """
    MATCH (g:Gene)
    WHERE
        toLower(toString(coalesce(g.symbol, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(g.name, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(g.gene_id, '')))
            CONTAINS $q
    RETURN properties(g) AS gene
    LIMIT $limit
    """

    # Search diseases by common name and identifier properties.
    disease_query = """
    MATCH (d:Disease)
    WHERE
        toLower(toString(coalesce(d.name, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(d.disease_name, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(d.disease_id, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(d.id, '')))
            CONTAINS $q
    RETURN properties(d) AS disease
    LIMIT $limit
    """

    # Search connected gene-disease pairs.
    # The relationship type is returned dynamically.
    association_query = """
    MATCH (g:Gene)-[r]-(d:Disease)
    WHERE
        toLower(toString(coalesce(g.symbol, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(g.name, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(g.gene_id, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(d.name, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(d.disease_name, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(d.disease_id, '')))
            CONTAINS $q
        OR toLower(toString(coalesce(d.id, '')))
            CONTAINS $q

    RETURN
        properties(g) AS gene,
        type(r) AS relationship,
        properties(d) AS disease
    LIMIT $limit
    """

    genes = run_query(gene_query, parameters)
    diseases = run_query(disease_query, parameters)
    associations = run_query(association_query, parameters)

    return {
        "query": q,
        "genes": genes,
        "diseases": diseases,
        "associations": associations,
        "counts": {
            "genes": len(genes),
            "diseases": len(diseases),
            "associations": len(associations)
        }
    }
    
@app.get("/api/graph")
def get_graph(
    q: str = Query(..., min_length=1, max_length=200)
):
    query = q.strip()

    if not query:
        raise HTTPException(
            status_code=400,
            detail="Enter a gene or disease to visualize.",
        )

    cypher = """
    MATCH (n)
    WHERE any(
        key IN keys(n)
        WHERE toLower(toString(n[key]))
              CONTAINS toLower($query)
    )
    WITH n
    LIMIT 15

    OPTIONAL MATCH (n)-[r]-(other)

    RETURN n, r, other
    """

    nodes = {}
    edges = {}

    def make_node(node):
        node_id = node.element_id
        props = dict(node)
        labels = list(node.labels)

        # Prefer meaningful names for graph labels.
        display_label = next(
            (
                str(props[key])
                for key in (
                    "name",
                    "symbol",
                    "gene_symbol",
                    "disease_name",
                    "label",
                    "id",
                )
                if props.get(key) is not None
            ),
            labels[0] if labels else "Node",
        )

        return {
            "data": {
                "id": node_id,
                "label": display_label,
                "labels": labels,
                "properties": props,
            }
        }

    try:
        with driver.session() as session:
            records = session.run(cypher, query=query)

            for record in records:
                node = record["n"]
                other = record["other"]
                relationship = record["r"]

                if node is not None:
                    nodes[node.element_id] = make_node(node)

                if other is not None:
                    nodes[other.element_id] = make_node(other)

                if relationship is not None:
                    source = relationship.start_node.element_id
                    target = relationship.end_node.element_id
                    rel_type = relationship.type

                    edge_id = relationship.element_id

                    edges[edge_id] = {
                        "data": {
                            "id": edge_id,
                            "source": source,
                            "target": target,
                            "label": rel_type,
                        }
                    }

        return {
            "query": query,
            "node_count": len(nodes),
            "edge_count": len(edges),
            "elements": {
                "nodes": list(nodes.values()),
                "edges": list(edges.values()),
            },
        }

    except Neo4jError as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Could not build graph: {exc}",
        )
