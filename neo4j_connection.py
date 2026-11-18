
import os
from neo4j import GraphDatabase

URI = "neo4j://127.0.0.1:7687"
USERNAME = "neo4j"
PASSWORD = os.environ.get("NEO4J_PASSWORD")

if not PASSWORD:
    raise ValueError(
        "Please set the NEO4J_PASSWORD environment variable."
    )

with GraphDatabase.driver(
    URI,
    auth=(USERNAME, PASSWORD)
) as driver:
    driver.verify_connectivity()
    print("Connected to Neo4j successfully!")

    with driver.session() as session:
        result = session.run(
            "RETURN 'GeneDiseaseKG is connected!' AS message"
        )
        print(result.single()["message"])
