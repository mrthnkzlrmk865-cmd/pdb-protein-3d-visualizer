# PDB Protein 3D Structural Visualizer & Analyzer

An advanced structural bioinformatics pipeline written in Python that programmatically retrieves verified 3D macromolecular crystal structures and analyzes atomic trajectories.

## Features
- **Live API Integration:** Dynamically fetches experimental molecular structures directly from the RCSB Protein Data Bank (PDB).
- **Biopython Parsing:** Utilizes `Bio.PDB` structural hierarchy to parse and extract specific atomic coordinate datasets (Alpha Carbons - CA).
- **Stereoscopic 3D Trajectory:** Models the precise peptide folding backbone using spatial coordinate matrices in Ångström unit scales.

## Visual Output
Below is the simulated 3D spatial trajectory of the Human Tumor Suppressor P53 tetramerization domain (PDB ID: 1OLG):

![Protein 3D Plot](p53_3d_backbone.png)

## Requirements
To execute this structural analyzer pipeline locally, install the core dependencies:
```bash
pip install biopython matplotlib numpy
