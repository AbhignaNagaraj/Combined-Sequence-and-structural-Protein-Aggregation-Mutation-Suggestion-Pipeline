**Protein Aggregation & Mutation Suggestion Pipeline**


**Overview**

This repository provides an automated protein developability analysis pipeline that integrates sequence-based and structure-based aggregation propensity analysis to identify high-confidence aggregation hotspots and generate rational solubility-enhancing mutation suggestions.

The pipeline combines intrinsic hydrophobic clustering from the amino-acid sequence with solvent accessibility analysis from a predicted 3D structure (AlphaFold or homology model), enabling informed mutation design for recombinant protein expression.


**Key Features**

Sequence-based aggregation detection using hydrophobic clustering

Structure-based aggregation detection using solvent-accessible surface area (SASA)

Integration of sequence and structural signals to identify high-confidence aggregation hotspots

Automated, biophysically sound mutation suggestions

Single unified CSV report per protein

Single PNG visualization summarizing aggregation risk

Compatible with AlphaFold and SWISS-MODEL PDB files

Uses only free, open-source Python libraries

**Scientific Rationale**
Sequence-based aggregation analysis identifies intrinsically aggregation-prone regions independent of folding, while structure-based analysis identifies surface-exposed hydrophobic residues capable of mediating aggregation in the folded protein.
By intersecting both signals, the pipeline reduces false positives and highlights actionable aggregation drivers suitable for rational mutagenesis.


**Input Requirements**

**1. Protein Sequence (FASTA)**
Single-sequence FASTA file:
input_protein.fasta


**2. Protein Structure (PDB)**
Predicted 3D structure from:
**AlphaFold Protein **
Structure Database

**SWISS-MODEL**
Example:
input_protein.pdb


**Output Files**

All results are written to the results/ directory.

**results/**

├── <ProteinID>_aggregation_mutation_analysis.csv

└── <ProteinID>_aggregation_mutation_analysis.png

CSV Output Columns

Column	Description

Position	Residue index

Residue	Wild-type amino acid

Sequence_Aggregation	Hydrophobic clustering flag

Structural_Aggregation	Surface exposure (SASA-based)

High_Confidence_Aggregation	Intersection of both methods


**Suggested_Mutation	Automated solubility-enhancing mutation**

**Mutation Strategy**

Mutation suggestions target surface-exposed hydrophobic residues using conservative substitutions designed to improve solubility while preserving structural integrity:

Aliphatic hydrophobics → Aspartic acid (charge + hydration)

Aromatic residues → Asparagine (polar, minimally disruptive)

Alanine → Serine (mild polarity increase)

This strategy is widely used in protein engineering and developability optimization.

**Installation**

**Requirements**

Python ≥ 3.8

**Dependencies**

pip install biopython pandas matplotlib freesasa

**Usage**

Run the pipeline from the project directory:

python3 Aggregation_Mutation.py

**Ensure the following files are present in the same directory:**

input_protein.fasta

input_protein.pdb


**Visualization**

The generated PNG plot highlights high-confidence aggregation residues (sequence + structure) along the protein sequence, providing a rapid visual overview of aggregation risk and mutation targets.


**Applications**

Recombinant protein expression optimization

Solubility engineering

Protein developability assessment

Pre-experimental mutation prioritization

Educational and training purposes in computational biology


**Limitations**

Structural analysis depends on the accuracy of the predicted 3D model

Sequence analysis assumes linear hydrophobic clustering

Experimental validation is required to confirm predicted improvements


**Author**

Developed as part of a technical bioinformatics assessment focused on protein aggregation analysis, structural evaluation, and rational mutation design.


**License**
This project is provided for academic, research, and evaluation purposes.

