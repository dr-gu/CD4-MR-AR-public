# CD4-MR-AR-public

Analysis code accompanying the CD4-MR-AR study: genetic prioritization
of CD4+ T cell-expressed candidate genes for allergic rhinitis, spanning
transcriptome-wide Mendelian randomization (MR) across 46 dynamic CD4+
T cell eQTL profiles, colocalization and SMR/HEIDI refinement,
cross-database replication, metabolite mediation screening, single-cell
and paired TCR contextualization, and structure-guided screening of
SLC25A46 and TLR1.

## Repository contents

| Path | Description |
|---|---|
| `scripts/` | One analysis script per Results subsection of the manuscript |
| `ENVIRONMENT.md` | Machine, R/Python versions and key packages used |
| `LICENSE` | MIT license |

### Scripts

| Script | Results | Analysis |
|---|---|---|
| `01_primary_mr.R` | 2.1; Figs. 1–2 | Soskic CD4+ dynamic eQTL instrument extraction, rsID mapping, harmonization with the AR GWAS (GCST90468131), MHC-region exclusion, LD clumping, Steiger filtering, primary MR (Wald/IVW/MR-RAPS), genomic landscape, volcano and forest plots |
| `02_coloc_smr_heidi.R` | 2.2; Fig. 3 | Bayesian colocalization (coloc), timepoint-resolved support, SMR + HEIDI heterogeneity test, top-25 bar chart, heatmap and MR-vs-SMR comparison |
| `03_crossdb_overlap.R` | 2.3; Fig. 4, Table 1 | Reference eQTL extraction, DICE and eQTLGen replication MR, three-database overlap (nine-gene set), direction concordance, timepoint-dependent trajectories |
| `04_metabolite_mediation.R` | 2.4; Fig. 5 | Gene-to-metabolite MR across 314 metabolite GWAS, metabolite-to-AR MR, indirect-effect mediation chains, metabolite coloc |
| `05_scrnaseq_tcr.R` | 2.5; Fig. 6, Figs. S1–S10 | Seurat single-cell workflows for the four datasets, cell annotation, MR-gene expression landscape, Ro/e enrichment and Fisher tests, pseudobulk reanalysis, miloR neighborhood DA, paired TCR clonotype analysis, evidence matrix |
| `06_structure_screening.R` | 2.6; Fig. 7, Table S13 | Representative excerpt of the GraphBAN virtual screening and ADMET-AI filtering pipeline, AutoDock Vina docking, GROMACS short MD, Figure 7 assembly |
| `R_functions_roe_analysis.R` | 2.5 helpers | Shared Ro/e, Fisher, miloR and figure-theme helpers sourced by `05_scrnaseq_tcr.R` |

## Data availability

The scripts read input data that are obtained from public repositories as
described in the manuscript Methods: the AR GWAS summary statistics
(GCST90468131), Soskic CD4+ T cell dynamic eQTL data, DICE and eQTLGen
eQTL resources, 314 metabolite GWAS from the GWAS Catalog, and the
single-cell RNA-seq datasets GSE200107 and GSE273975 from GEO.
Intermediate processed files are not redistributed
with this repository; place the downloaded inputs under `data/` as
indicated at the top of each script.

## Citation

Please cite the published article.
