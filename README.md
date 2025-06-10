**Molecular Representations Generator**

This repository contains four Python scripts designed to generate six types of molecular representations from SMILES strings. The generated representations can be used for cheminformatics, QSAR modeling, or other machine learning tasks in computational chemistry.

Molecular Representations
The following six molecular representations are supported:

Fingerprints:
- Morgan fingerprints

- MACCS keys

- RDKit fingerprints

Descriptors:
- Mordred descriptors

- RDKit descriptors

- Selected 13 RDKit descriptors (a curated subset of common physicochemical properties)

**Functionality**
Each script processes molecular structures provided in SMILES format, extracts the specified representation, and automatically generates a merged output table that includes:

The computed molecular representation

The associated experimental parameters from the input dataset (e.g., biological activity)

This streamlined workflow ensures that both molecular features and experimental data are ready for downstream analysis or machine learning.

**Input**
The input is a CSV file containing at least one column with SMILES strings and additional columns with experimental data.

**Output**
Each script saves a cleaned and processed dataset, ready for further modeling or statistical analysis.
