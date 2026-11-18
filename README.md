# 🧬 Gene–Disease Knowledge Graph
### A Big Data Analytics and Network-Based Approach to Understanding Gene–Disease Associations

## 📌 Project Overview

The **Gene–Disease Knowledge Graph** is a Big Data Analytics project designed to explore and represent relationships between human genes and diseases using publicly available biological datasets.

The project integrates large-scale biomedical data from **Open Targets** and the **GWAS Catalog** to construct a structured gene–disease association network. By leveraging Python-based data processing and Neo4j graph database technology, the project aims to organize complex biological relationships into an interconnected graph representation.

Knowledge graphs provide an effective way to represent biological entities and their relationships, making it easier to explore gene–disease associations and identify potential connections within biomedical datasets.

## 🎯 Objectives

- Integrate gene-related and disease-association information from multiple biomedical datasets.
- Process and transform large-scale biological datasets using Python.
- Explore gene–disease relationships from Open Targets and GWAS Catalog.
- Represent biological entities and their associations using a graph database.
- Implement a Neo4j-based knowledge graph for storing and exploring biological relationships.
- Demonstrate the application of Big Data technologies in computational biology and biomedical research.

## 📊 Data Sources

The project utilizes publicly available biomedical resources:

### 1. Open Targets
Open Targets provides comprehensive information about relationships between genes, diseases, and therapeutic targets.

Website: https://platform.opentargets.org/

### 2. GWAS Catalog
The GWAS Catalog provides curated genome-wide association study data, including associations between genetic variants and different traits and diseases.

Website: https://www.ebi.ac.uk/gwas/

## 🛠️ Technologies Used

| Technology | Purpose |
|---|---|
| Python | Data processing and integration |
| Pandas | Data manipulation and preprocessing |
| PyArrow | Handling large Parquet datasets |
| Neo4j | Graph database management |
| Cypher | Graph querying and relationship exploration |
| CSV | Structured data exchange |
| Parquet | Efficient storage of large-scale datasets |

## ⚙️ Project Workflow

The project follows these major stages:

**1. Data Collection**
- Obtaining gene and disease-related information from Open Targets.
- Collecting genetic association information from the GWAS Catalog.

**2. Data Preprocessing**
- Loading Parquet and TSV datasets using Python.
- Cleaning and transforming relevant biological information.
- Extracting gene identifiers, disease information, and association records.

**3. Data Integration**
- Combining relevant information from multiple biological data sources.
- Organizing the extracted information into structured datasets.
- Preparing the data for graph database ingestion.

**4. Knowledge Graph Construction**
- Creating graph entities representing genes and diseases.
- Establishing relationships based on available gene–disease association records.
- Importing structured data into Neo4j.

**5. Graph Exploration**
- Using Cypher queries to explore relationships.
- Investigating connections between genes and associated diseases.
- Representing biological associations through a graph-based structure.

## 🏗️ Knowledge Graph Architecture

The conceptual structure of the knowledge graph consists of:

**Nodes:**
- Gene
- Disease
- Genetic Association

**Relationships:**
- Gene → ASSOCIATED_WITH → Disease
- Genetic Association → SUPPORTS → Gene–Disease Association

*The final node labels and relationship types depend on the implemented graph schema.*

## 📁 Project Structure

```text
Gene-Disease-Knowledge-Graph/
│
├── data/
│   ├── genes.csv
│   ├── opentargets/
│   └── gwas_catalog/
│
├── scripts/
│   ├── data_preprocessing.py
│   ├── data_integration.py
│   └── neo4j.py
│
├── notebooks/
│   └── exploratory_analysis.ipynb
│
├── README.md
└── requirements.txt
```

*The directory structure is illustrative and should be adjusted to match the actual repository.*

## 💻 Installation and Setup

### 1. Clone the Repository

```bash
git clone https://github.com/Uzmahaseeb4/Gene-Disease-Knowledge-Graph.git
```

### 2. Navigate to the Project Directory

```bash
cd Gene-Disease-Knowledge-Graph
```

### 3. Create a Virtual Environment

```bash
python -m venv env
```

Activate the environment:

Windows:
```bash
env\Scripts\activate
```

### 4. Install Dependencies

```bash
pip install pandas pyarrow neo4j
```

### 5. Configure Neo4j

- Install and launch Neo4j Desktop or Neo4j Community Edition.
- Create a database.
- Configure the database connection credentials.
- Connect the Python application to the Neo4j database.

## 🔬 Applications

This project demonstrates how graph-based data representation can support:

- Gene–disease association exploration.
- Biomedical knowledge representation.
- Integration of heterogeneous biological datasets.
- Computational biology research.
- Network-based analysis of biological relationships.
- Future development of disease-gene discovery and prioritization systems.

## 🚀 Future Scope

Potential extensions include:

- Integration of additional biomedical databases such as DisGeNET and UniProt.
- Incorporation of protein–protein interaction networks.
- Development of graph-based disease association analysis.
- Identification of shared genes across multiple diseases.
- Application of graph algorithms and graph machine learning.
- Development of an interactive knowledge graph visualization interface.

## 👩‍💻 Author

**Uzma Haseeb**

M.Sc. Data Science and Bio AI  
Digital University Kerala

- GitHub: https://github.com/Uzmahaseeb4
- LinkedIn: https://www.linkedin.com/in/uzma-ba5195318/

## 📚 Acknowledgements

We acknowledge the Open Targets Platform, GWAS Catalog, Neo4j, and the open-source scientific community for providing the resources and technologies that support this project.

## 📜 License

This project is intended for academic and educational purposes.
