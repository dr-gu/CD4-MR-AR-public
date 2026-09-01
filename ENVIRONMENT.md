# Compute Environment

All analyses were run on a single workstation:

- **Hardware**: Apple M2 Pro, 32 GB unified memory, arm64
- **OS**: macOS 26.6.2
- **Package management**: Miniconda (Python 3.11 envs), Homebrew

## Software and versions

| Software | Version | Role |
|---|---|---|
| R | 4.5.3 | MR and single-cell pipelines |
| data.table | 1.18.4 | large-table I/O |
| TwoSampleMR | 0.6.29 | MR analysis |
| coloc | 5.2.3 | colocalization |
| mr.raps | 0.4.3 | MR-RAPS estimator |
| arrow | 23.0.1.2 | parquet eQTL input |
| ggplot2 / ggrepel | 4.0.3 / 0.9.8 | visualization |
| dplyr / tidyr / patchwork | 1.2.1 / 1.3.2 / 1.3.2 | data handling and figure assembly |
| VennDiagram | 1.8.2 | gene overlap |
| Seurat | 5.5.0 | single-cell RNA-seq analysis |
| harmony | 1.2.4 | single-cell batch integration |
| sctransform | 0.4.3 | single-cell normalization |
| celldex / SingleR | 1.20.0 / 2.12.0 | reference-based cell annotation |
| Python | 3.11.15 | docking, screening, and figure scripts |
| RDKit | 2024.03.6 | ligand handling |
| meeko | 0.7.1 | docking input preparation |
| matplotlib | 3.10.9 | plotting |
| numpy / pandas | 2.4.6 / 2.3.3 | numerical and tabular data |
| AutoDock Vina | 1.2.6 (conda-forge) | molecular docking |
| PyTorch | 2.3.1 | GraphBAN screening models |
| ADMET-AI | 2.0.1 | ADMET filtering |
| GROMACS | 2025.4 (conda-forge) | molecular dynamics |
| PyMOL | 3.x (open-source, conda-forge) | structure visualization |
| pandoc | 3.10 | document conversion |
| XeTeX | TeX Live 2026 | PDF typesetting |
