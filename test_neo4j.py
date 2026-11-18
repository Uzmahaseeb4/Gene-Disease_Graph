
import os
from neo4j import GraphDatabase

URI = "neo4j://127.0.0.1:7687"
USERNAME = "neo4j"
PASSWORD = os.getenv("NEO4J_PASSWORD")

if not PASSWORD:
    raise RuntimeError(
        "Set the NEO4J_PASSWORD environment variable first."
    )

def main():
    with GraphDatabase.driver(
        URI, auth=(USERNAME, PASSWORD)
    ) as driver:
        driver.verify_connectivity()
        print("Connection successful!")

        with driver.session() as session:
            result = session.run(
                "RETURN 1 AS test"
            )
            print("Query result:", result.single()["test"])

if __name__ == "__main__":
    main()
