# =============================================================================
# 06_structure_screening.R
# Results 2.6: Focused computational follow-up prioritized representative
# binders for SLC25A46 and TLR1 (paper Figure 7, Supplementary Table S13)
# Paper: Genetic prioritization of CD4+ T cell-expressed candidate genes for
#        allergic rhinitis (Biomedicines)
# Merged and sanitized from: drug_config.R, drug_01-drug_10 analysis scripts
# =============================================================================
# Notes:
#  - Run from the repository root; each script assumes cwd == repo root.
#  - External engines (GraphBAN, ADMET-AI, AutoDock Vina, GROMACS, PyMOL)
#    run in separate conda environments; binary paths are placeholders.
#  - The screening library (ZINC-Pin1) is obtained from public sources.
#  - Repetitive blocks are elided and marked "# (elided: ...)".
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Configuration (all binary paths are placeholders)
# -----------------------------------------------------------------------------
# Two Python environments are required because DGL (GraphBAN) and admet_ai
# have incompatible PyTorch requirements:
PYTHON_BIN       <- "<path-to-conda-env>/bin/python"  # graphban-clean: torch 2.3.1 + DGL + ESM + transformers
PYTHON_ADMET_BIN <- "<path-to-conda-env>/bin/python"  # admet-ai: torch 2.8+ + admet_ai
GMX_BIN          <- "<path-to-md-env>/bin/gmx"        # md-clean: GROMACS
# GraphBAN workflow scripts (screening, ADMET filtering, receptor/box
# preparation, Vina batch wrapper) and the public ZINC-Pin1 library:
GRAPHBAN_AUTO <- "<path-to-graphban-repo>/automation"
ZINC_PIN1     <- "data/drug/ZINC-Pin1.csv"

# Repo-relative output locations
DRUG_RESULTS <- "results/drug"
DRUG_FIGURES <- "figures/drug"
DRUG_TABLES  <- "tables/drug"

# Per-target output directories
create_drug_dirs <- function() {
  sub <- c("predictions/biosnap", "predictions/kiba", "admet",
           "docking/receptor", "docking/receptor_pdbqt", "docking/box",
           "docking/ligands", "docking/ligands_pdbqt", "docking/homolog", "md")
  for (target in c("TLR1", "SLC25A46"))
    for (s in sub)
      dir.create(file.path(DRUG_RESULTS, s, target), recursive = TRUE, showWarnings = FALSE)
  invisible(lapply(c(DRUG_FIGURES, DRUG_TABLES), dir.create, recursive = TRUE, showWarnings = FALSE))
}
create_drug_dirs()

# --- AutoDock Vina parameters ---
VINA_EXHAUSTIVENESS <- 32
VINA_NUM_MODES      <- 9
VINA_ENERGY_RANGE   <- 4

# --- GROMACS MD parameters (short explicit-solvent production) ---
MD_NSTEPS_EM   <- 50000
MD_NSTEPS_NVT  <- 50000
MD_NSTEPS_NPT  <- 50000
MD_NSTEPS_PROD <- 500000      # 500000 x 0.002 ps = 1 ns production
MD_DT          <- 0.002
MD_TEMP        <- 298.15
MD_PRESSURE    <- 1.0
MD_FF          <- "amber99sb-ildn"
MD_WATER       <- "tip3p"
MD_SALT        <- "0.15"

library(data.table)
library(ggplot2)
library(scales)
library(circlize)
library(fmsb)
library(patchwork)
library(jsonlite)
library(magick)
library(png)

# -----------------------------------------------------------------------------
# 2. Virtual screening (GraphBAN) and pharmacokinetic filtering (ADMET-AI)
# -----------------------------------------------------------------------------
# Screening library: ZINC-Pin1, 498,910 compounds. SMILES deduplication leaves
# 249,455 unique compounds; 237,326 pass Lipinski + PAINS-like drug-likeness
# filtering and enter scoring.
# Compound features: ChemBERTa molecular embeddings; protein features: ESM-1b
# embeddings of UniProt sequences (Homo sapiens, taxon 9606).
# Four GraphBAN models are scored (BioSNAP interaction + KIBA affinity);
# the consensus score is the mean across the four models.
genes <- c("TLR1", "SLC25A46")

# --- GraphBAN screening (graphban-clean environment) ---
cmd <- paste(shQuote(PYTHON_BIN),
  file.path(GRAPHBAN_AUTO, "graphban_target_screen.py"),
  "--genes", paste(genes, collapse = " "),
  "--organism-id", "9606",
  "--zinc-source", shQuote(ZINC_PIN1),
  "--output-dir", shQuote(DRUG_RESULTS),
  "--chunk-size", "10000",
  "--datasets", "biosnap", "kiba",
  "--device", "auto")
message("Running GraphBAN screening for: ", paste(genes, collapse = ", "))
system(cmd)

# --- ADMET-AI filtering (admet-ai environment) ---
# Nine ADMET endpoints: mutagenicity, carcinogenicity, clinical toxicity,
# drug-induced liver injury, hERG inhibition, oral bioavailability, half-life,
# plasma protein binding, logP. 164 compounds per target pass all nine rules.
for (gene in genes) {
  cmd <- paste(shQuote(PYTHON_ADMET_BIN),
    file.path(GRAPHBAN_AUTO, "run_admet_filter.py"),
    "--gene", gene,
    "--predictions-root", file.path(DRUG_RESULTS, "predictions"),
    "--output-dir", file.path(DRUG_RESULTS, "admet", gene),
    "--batch-size", "1000")
  message("Running ADMET-AI filtering for: ", gene)
  system(cmd)
}

# -----------------------------------------------------------------------------
# 3. Molecular docking (AutoDock Vina)
# -----------------------------------------------------------------------------
# Four docking lines: TLR1 LRR, TLR1 TIR, SLC25A46 full-length AlphaFold
# pocket (Strategy A), SLC25A46 mitochondrial carrier family domain
# (Strategy B). Receptor models are AlphaFold structures (UniProt Q15399 for
# TLR1, Q96AG3 for SLC25A46). The TLR1 TIR-domain pocket outperformed the LRR
# region (lead ZINC283870: -5.856 kcal/mol vs -4.147 kcal/mol); the SLC25A46
# lead is ZINC36371100. One representative line (TLR1 TIR) is shown in full.

run_cmd <- function(cmd, desc) {
  message("Running: ", desc)
  if (system(cmd, intern = FALSE) != 0) stop("Command failed: ", desc)
}

# --- Download AlphaFold models (cached locally) ---
download_alphafold <- function(uniprot, output_path) {
  if (file.exists(output_path)) return()
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot, "-F1-model_v6.pdb")
  download.file(url, output_path, method = "curl")
}

# --- Trim receptor to the docking domain ---
# Residue ranges follow the domain annotations in the manuscript methods;
# update them for the local AlphaFold models.
trim_receptor <- function(input_pdb, output_pdb, start_res, end_res,
                           metadata_json, gene, source_url, notes) {
  cmd <- paste(shQuote(PYTHON_BIN),
    file.path(GRAPHBAN_AUTO, "prepare_receptor_from_alphafold.py"),
    "--input-pdb", shQuote(input_pdb), "--output-pdb", shQuote(output_pdb),
    "--start", start_res, "--end", end_res,
    "--metadata-json", shQuote(metadata_json), "--gene", gene,
    "--source-url", shQuote(source_url), "--notes", shQuote(notes))
  run_cmd(cmd, paste("Trim receptor:", gene, start_res, "-", end_res))
}

# Representative call: TLR1 TIR domain
download_alphafold("Q15399", "data/drug/alphafold/AF-Q15399-F1-model_v6.pdb")
trim_receptor(
  input_pdb     = "data/drug/alphafold/AF-Q15399-F1-model_v6.pdb",
  output_pdb    = file.path(DRUG_RESULTS, "docking/TLR1/receptor/TLR1_TIR.pdb"),
  start_res     = "<domain-start-residue>", end_res = "<domain-end-residue>",
  metadata_json = file.path(DRUG_RESULTS, "docking/TLR1/receptor/TLR1_TIR_metadata.json"),
  gene = "TLR1", source_url = "<alphafold-entry-url>", notes = "TIR domain pocket")
# (elided: TLR1 LRR trim and the two SLC25A46 receptor trims)

# --- Derive docking box from a homologous complex ---
derive_box <- function(target_pdb, homolog_pdb, homolog_chain, homolog_ligand,
                        output_json, output_config, receptor_pdbqt,
                        padding = 8.0, min_size = 20.0, site_dist = 6.0) {
  cmd <- paste(shQuote(PYTHON_BIN),
    file.path(GRAPHBAN_AUTO, "derive_vina_box_from_homolog.py"),
    "--target-pdb", shQuote(target_pdb), "--homolog-pdb", shQuote(homolog_pdb),
    "--homolog-chain", homolog_chain, "--homolog-ligand", homolog_ligand,
    "--site-distance", site_dist, "--padding", padding, "--minimum-size", min_size,
    "--output-json", shQuote(output_json), "--output-config", shQuote(output_config),
    "--receptor-pdbqt", shQuote(receptor_pdbqt))
  run_cmd(cmd, paste("Derive box:", basename(target_pdb), "<-", basename(homolog_pdb)))
}

# Representative call: TLR1 TIR pocket box
derive_box(
  target_pdb   = file.path(DRUG_RESULTS, "docking/TLR1/receptor/TLR1_TIR.pdb"),
  homolog_pdb  = "data/drug/homolog/<homolog-complex>.pdb",
  homolog_chain = "<chain-id>", homolog_ligand = "<ligand-name>",
  output_json  = file.path(DRUG_RESULTS, "docking/TLR1/box/TLR1_TIR_box.json"),
  output_config = file.path(DRUG_RESULTS, "docking/TLR1/box/TLR1_TIR_box.cfg"),
  receptor_pdbqt = file.path(DRUG_RESULTS, "docking/TLR1/receptor_pdbqt/TLR1_TIR.pdbqt"))
# (elided: box derivation for TLR1 LRR and the two SLC25A46 strategies)

# --- Prepare ligands from the ADMET-passed set (top 20 per target) ---
prepare_ligands <- function(gene) {
  cmd <- paste(shQuote(PYTHON_BIN),
    file.path(GRAPHBAN_AUTO, "prepare_ligands_from_admet.py"),
    "--input-csv", shQuote(file.path(DRUG_RESULTS, "admet", gene, paste0(gene, "_admet_passed.csv"))),
    "--output-dir", shQuote(file.path(DRUG_RESULTS, "docking", gene, "ligands")),
    "--top-k", "20", "--ranking-column", "rank_after_admet")
  run_cmd(cmd, paste("Prepare ligands:", gene))
}

# --- Convert SDF to PDBQT with meeko (ligands and receptors) ---
convert_ligands_pdbqt <- function(gene) {
  ligand_dir <- file.path(DRUG_RESULTS, "docking", gene, "ligands")
  pdbqt_dir <- file.path(DRUG_RESULTS, "docking", gene, "ligands_pdbqt")
  dir.create(pdbqt_dir, recursive = TRUE, showWarnings = FALSE)
  for (sdf in list.files(file.path(ligand_dir, "sdf"), pattern = "\\.sdf$", full.names = TRUE)) {
    base <- tools::file_path_sans_ext(basename(sdf))
    pdbqt_out <- file.path(pdbqt_dir, paste0(base, ".pdbqt"))
    if (file.exists(pdbqt_out)) next
    py_cmd <- paste(shQuote(PYTHON_BIN), "-c", shQuote(paste0(
"from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
mol = Chem.SDMolSupplier('", sdf, "', removeHs=False)[0]
prep = MoleculePreparation()
mol_setup = prep.prepare(mol)[0]
pdbqt_string, is_ok = PDBQTWriterLegacy.write_string(mol_setup)
if is_ok:
    with open('", pdbqt_out, "', 'w') as f:
        f.write(pdbqt_string)
    print('OK ", base, "')
else:
    print('FAILED ", base, "')
")))
    system(py_cmd, intern = FALSE)
  }
}
# (elided: analogous meeko conversion of receptor PDB structures to PDBQT)

prepare_ligands("TLR1")
convert_ligands_pdbqt("TLR1")
# (elided: ligand preparation and PDBQT conversion for SLC25A46)

# --- Run Vina docking per line ---
run_vina_line <- function(receptor_pdbqt, ligand_dir, output_dir,
                           center_x, center_y, center_z,
                           size_x = 20, size_y = 20, size_z = 20) {
  cmd <- paste(shQuote(PYTHON_BIN),
    file.path(GRAPHBAN_AUTO, "run_vina_batch.py"),
    "--receptor", shQuote(receptor_pdbqt), "--ligand-dir", shQuote(ligand_dir),
    "--output-dir", shQuote(output_dir),
    "--center-x", center_x, "--center-y", center_y, "--center-z", center_z,
    "--size-x", size_x, "--size-y", size_y, "--size-z", size_z,
    "--exhaustiveness", VINA_EXHAUSTIVENESS,
    "--num-modes", VINA_NUM_MODES, "--energy-range", VINA_ENERGY_RANGE)
  run_cmd(cmd, paste("Vina docking:", basename(output_dir)))
}

# Representative call: TLR1 TIR pocket; box center from the derived box config
run_vina_line(
  receptor_pdbqt = file.path(DRUG_RESULTS, "docking/TLR1/receptor_pdbqt/TLR1_TIR.pdbqt"),
  ligand_dir     = file.path(DRUG_RESULTS, "docking/TLR1/ligands_pdbqt"),
  output_dir     = file.path(DRUG_RESULTS, "docking/TLR1/vina_top20_TIR"),
  center_x = "<box-center-x>", center_y = "<box-center-y>", center_z = "<box-center-z>",
  size_x = 20, size_y = 20, size_z = 20)
# (elided: Vina runs for TLR1 LRR, SLC25A46 full-length and MCF strategies)

# -----------------------------------------------------------------------------
# 4. Molecular dynamics (GROMACS, 1 ns explicit solvent)
# -----------------------------------------------------------------------------
# Four MD systems follow the docking lines. Each system is energy-minimized,
# equilibrated (NVT, NPT), then run for 500000 x 0.002 ps = 1 ns production in
# explicit solvent (amber99sb-ildn, TIP3P, 0.15 M salt, 298.15 K), with three
# replicate simulations per system. Reported values for the TLR1 TIR
# replicates: mean backbone RMSD 0.202 nm, mean 296 protein-ligand contacts
# per frame.
MD_SYSTEMS <- c("TLR1_LRR_rank001", "TLR1_TIR_rank001",
                "SLC25A46_full_rank001", "SLC25A46_MCF_rank001")

run_cmd_md <- function(cmd, desc) {
  message("Running: ", desc)
  if (system(cmd, intern = FALSE) != 0) warning("Command returned non-zero: ", desc)
}

# --- Lead selection: docking rank 001 of each line ---
select_leads <- function() {
  docking_lines <- c("TLR1/vina_top20_LRR", "TLR1/vina_top20_TIR",
                     "SLC25A46/vina_top20_full", "SLC25A46/vina_top20_MCF")
  leads <- list()
  for (line in docking_lines) {
    summary_file <- file.path(DRUG_RESULTS, "docking", line, "docking_summary.csv")
    if (!file.exists(summary_file)) { warning("Docking summary not found: ", summary_file); next }
    dt <- fread(summary_file)
    setorder(dt, best_affinity_kcal_mol)
    lead <- dt[1]
    leads[[line]] <- lead
    message(sprintf("%s: lead = %s, affinity = %.3f kcal/mol",
                    line, lead$ligand_id, lead$best_affinity_kcal_mol))
  }
  leads
}

# --- GROMACS analysis commands (one representative system; run per system) ---
run_basic_md_analysis <- function(sys_dir) {
  gmx <- GMX_BIN
  analysis_dir <- file.path(sys_dir, "analysis_basic")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  tpr <- file.path(sys_dir, "md.tpr"); xtc <- file.path(sys_dir, "md.xtc")
  if (!file.exists(tpr) || !file.exists(xtc)) { warning("Missing tpr/xtc for: ", sys_dir); return() }
  # Backbone RMSD
  xvg <- file.path(analysis_dir, "backbone_rmsd.xvg")
  if (!file.exists(xvg))
    run_cmd_md(sprintf('echo "Backbone Backbone" | %s rms -s %s -f %s -fit rot+trans -o %s -tu ns',
                       gmx, tpr, xtc, xvg), paste("RMSD:", basename(sys_dir)))
  # C-alpha RMSF
  xvg <- file.path(analysis_dir, "ca_rmsf_residue.xvg")
  if (!file.exists(xvg))
    run_cmd_md(sprintf('echo "C-alpha C-alpha" | %s rmsf -s %s -f %s -res -o %s',
                       gmx, tpr, xtc, xvg), paste("RMSF:", basename(sys_dir)))
  # Radius of gyration
  xvg <- file.path(analysis_dir, "protein_rg.xvg")
  if (!file.exists(xvg))
    run_cmd_md(sprintf('echo "Protein" | %s gyrate -s %s -f %s -o %s',
                       gmx, tpr, xtc, xvg), paste("Rg:", basename(sys_dir)))
  # SASA
  xvg <- file.path(analysis_dir, "protein_sasa.xvg")
  if (!file.exists(xvg))
    run_cmd_md(sprintf('echo "Protein" | %s sasa -s %s -f %s -o %s',
                       gmx, tpr, xtc, xvg), paste("SASA:", basename(sys_dir)))
  # Protein-ligand minimum distance
  xvg <- file.path(analysis_dir, "protein_ligand_mindist.xvg")
  if (!file.exists(xvg))
    run_cmd_md(sprintf('echo "Protein Ligand" | %s mindist -s %s -f %s -od %s',
                       gmx, tpr, xtc, xvg), paste("PL distance:", basename(sys_dir)))
  # Protein-ligand contacts (0.45 nm cutoff)
  xvg <- file.path(analysis_dir, "protein_ligand_contacts.xvg")
  if (!file.exists(xvg))
    run_cmd_md(sprintf('echo "Protein Ligand" | %s mindist -s %s -f %s -on %s -d 0.45',
                       gmx, tpr, xtc, xvg), paste("PL contacts:", basename(sys_dir)))
}

# --- Hydrogen-bond analysis (MDAnalysis) ---
run_hbond_analysis <- function(sys_dir) {
  cmd <- paste(shQuote(PYTHON_BIN),
    file.path(GRAPHBAN_AUTO, "full_run_hybrid", "md", "hbond_analysis_mdanalysis.py"),
    "--topology", file.path(sys_dir, "md.tpr"),
    "--trajectory", file.path(sys_dir, "md.xtc"),
    "--output-dir", file.path(sys_dir, "analysis_basic"))
  run_cmd_md(cmd, paste("H-bond analysis:", basename(sys_dir)))
}

# --- Run all analyses on all systems ---
run_all_md_analyses <- function() {
  for (sys in MD_SYSTEMS) {
    sys_dir <- file.path(DRUG_RESULTS, "md", sys)
    if (!dir.exists(sys_dir)) { warning("System directory not found: ", sys_dir); next }
    run_basic_md_analysis(sys_dir)
    run_hbond_analysis(sys_dir)
  }
}

# -----------------------------------------------------------------------------
# 5. Representative visualization panels
# -----------------------------------------------------------------------------
# The full supplementary panel set (screening cascade, UMAP, GraphBAN
# model-consistency scatter, ADMET ridge distributions, ADMET rule-failure
# plot, GraphBAN-vs-Vina convergence, RMSF ridge, MD radar, candidate heatmap,
# Tanimoto network, MD metrics table) is elided; one panel per pattern is kept
# below. All panels are saved as PDF under figures/drug/.

theme_nature <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.3),
        axis.ticks = element_line(color = "black", linewidth = 0.3),
        legend.position = "bottom")

# --- Panel: dual-domain top-10 docking comparison (TLR1 TIR vs LRR shown) ---
generate_viz_binding_compare <- function(target, domains) {
  all_data <- rbindlist(lapply(names(domains), function(d) {
    f <- file.path(DRUG_RESULTS, "docking", target, domains[[d]], "docking_summary.csv")
    if (!file.exists(f)) return(NULL)
    dt <- fread(f); dt[, domain := d]; dt
  }))
  if (nrow(all_data) == 0) return()
  all_data[, rank := frank(best_affinity_kcal_mol), by = domain]
  top10 <- all_data[rank <= 10]
  p <- ggplot(top10, aes(x = reorder(ligand_id, -best_affinity_kcal_mol),
                          y = best_affinity_kcal_mol, fill = domain)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("#C65D3B", "#2F6B5F", "#E8926B", "#5BAF9A")) +
    coord_flip() + theme_nature +
    labs(x = "Ligand", y = "Vina best affinity (kcal/mol)", fill = "Domain")
  ggsave(file.path(DRUG_FIGURES, paste0("viz_", target, "_binding_compare.pdf")),
         p, width = 10, height = 6)
}
generate_viz_binding_compare("TLR1", c("LRR domain" = "vina_top20_LRR",
                                       "TIR domain" = "vina_top20_TIR"))
# (elided: SLC25A46 full vs MCF binding compare call)

# --- Panel: per-residue contact ring (circlize) ---
# contact_data: data.table(residue, contact_type, count)
# contact_type: H_bond, hydrophobic, pi_stack, salt_bridge
generate_viz_contact_circle <- function(contact_data, target, domain_label) {
  if (missing(contact_data) || nrow(contact_data) == 0) return()
  safe_label <- gsub("[ /]", "_", domain_label)
  pdf(file.path(DRUG_FIGURES,
      paste0("viz_", target, "_", safe_label, "_contact_circle.pdf")),
      width = 8, height = 8)
  chordDiagram(contact_data, annotationTrack = c("grid", "name"),
    grid.col = c(H_bond = "#2F6B5F", hydrophobic = "#C65D3B",
                 pi_stack = "#E8926B", salt_bridge = "#5BAF9A"),
    transparency = 0.3)
  title(paste(target, domain_label, "- residue-ligand contacts"))
  circos.clear()
  dev.off()
}

# --- Panel: MD stability traces (backbone RMSD time series) ---
read_xvg <- function(path) {
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[#@]", lines)]
  if (length(data_lines) == 0) return(NULL)
  fread(text = paste(data_lines, collapse = "\n"))
}
read_rmsd <- function(md_map) {
  rbindlist(lapply(names(md_map), function(nm) {
    f <- file.path(DRUG_RESULTS, "md", md_map[nm], "analysis_basic", "backbone_rmsd.xvg")
    if (!file.exists(f)) return(NULL)
    dt <- read_xvg(f)
    if (is.null(dt) || ncol(dt) < 2) return(NULL)
    setnames(dt, 1:2, c("time_ps", "rmsd_nm"))
    dt[, c("time_ns", "system") := list(time_ps / 1000, nm)]
    dt
  }))
}

MD_LEAD_MAP <- c("TLR1 TIR" = "TLR1_TIR_rank001",
                 "SLC25A46" = "SLC25A46_full_rank007")

generate_viz_md_rmsd <- function() {
  rmsd_all <- read_rmsd(MD_LEAD_MAP)
  if (nrow(rmsd_all) == 0) { message("No RMSD data"); return() }
  p <- ggplot(rmsd_all, aes(x = time_ns, y = rmsd_nm, color = system)) +
    geom_line(linewidth = 0.4) +
    scale_color_manual(values = c("TLR1 TIR" = "#2F6B5F", "SLC25A46" = "#C65D3B")) +
    theme_nature +
    labs(x = "Time (ns)", y = "Backbone RMSD (nm)", color = "System")
  ggsave(file.path(DRUG_FIGURES, "viz_md_rmsd_ts.pdf"), p, width = 10, height = 5)
}

# --- Panel: lead compound ADMET radar (normalized six-property profile) ---
generate_viz_admet_radar <- function() {
  lead_data <- list()
  for (target in c("TLR1", "SLC25A46")) {
    f <- file.path(DRUG_RESULTS, "admet", target, paste0(target, "_admet_passed.csv"))
    if (!file.exists(f)) next
    lead_data[[target]] <- fread(f)[1]
  }
  if (length(lead_data) < 1) { message("No lead data for ADMET radar"); return() }
  radar_mat <- do.call(rbind, lapply(names(lead_data), function(nm) {
    row <- lead_data[[nm]]
    col_val <- function(name) { v <- suppressWarnings(as.numeric(row[[name]])); if (!is.na(v)) v else NA_real_ }
    data.frame(MW = col_val("molecular_weight") / 500,
               logP = as.numeric(row$logP) / 5,
               TPSA = col_val("tpsa") / 140,
               Safety = 1 - as.numeric(row$AMES),
               Bioavail = as.numeric(row$Bioavailability_Ma),
               PPBR = as.numeric(row$PPBR_AZ) / 100,
               row.names = nm)
  }))
  radar_mat <- rbind(rep(1, 6), rep(0, 6), radar_mat)
  pdf(file.path(DRUG_FIGURES, "viz_admet_radar.pdf"), width = 7, height = 7)
  radarchart(radar_mat, axistype = 1,
    pcol = c("#C65D3B", "#2F6B5F"), plwd = 2, plty = 1,
    pfcol = alpha(c("#C65D3B", "#2F6B5F"), 0.2),
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    vlabels = c("MW", "logP", "TPSA", "Safety", "Bioavail", "PPBR"),
    vlcex = 0.8, title = "Lead compound ADMET radar")
  legend("topright", legend = names(lead_data), bty = "n",
         col = c("#C65D3B", "#2F6B5F"), lwd = 2)
  dev.off()
}

# -----------------------------------------------------------------------------
# 6. Paper Figure 7 assembly
# -----------------------------------------------------------------------------
# Figure 7 comprises 12 panels (a-l): SLC25A46 (a-f) and TLR1 (g-l), covering
# 3D binding-site pose, 2D ligand interaction diagram, domain docking
# comparison, top-5 docked poses, docking-to-MD overlay, and MD traces.
# Structure panels are PyMOL / RDKit renders; data panels are composed with
# patchwork. The composite layout is representative; per-panel duplication
# is elided. Output: figures/drug/Fig7.pdf (180 x 200 mm PDF).

# Helper: convert a square PNG (PyMOL render) or PDF (2D diagram) to an
# exact-size PDF panel
img_to_pdf_panel <- function(src, out_path, width_mm = 52.5, density = 600, is_pdf = FALSE) {
  img <- if (is_pdf) image_read_pdf(src, density = density) else image_read(src, density = density)
  img <- image_trim(img)
  side <- max(image_info(img)$width, image_info(img)$height)
  img <- image_extent(img, geometry = paste0(side, "x", side), gravity = "center", color = "white")
  target_px <- round(width_mm / 25.4 * density)
  img <- image_resize(img, paste0(target_px, "x", target_px))
  tmp <- tempfile(fileext = ".png")
  image_write(img, tmp, format = "png", density = paste0(density, "x", density))
  arr <- png::readPNG(tmp)
  unlink(tmp)
  pdf(out_path, width = width_mm / 25.4, height = width_mm / 25.4)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", bg = "white")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rasterImage(arr, 0, 0, 1, 1)
  box(col = "black", lwd = 0.5)
  dev.off()
}

# Representative render-embeds: TLR1 TIR 3D binding site (PyMOL) and 2D
# interaction diagram. SLC25A46 panels and remaining panels follow the pattern.
img_to_pdf_panel("figures/drug/viz_TLR1_TIR_binding_site.png",
                 "figures/drug/Fig7_panel_TLR1_TIR_3D.pdf")
img_to_pdf_panel("figures/drug/viz_TLR1_TIR_domain_02_ligand2d.pdf",
                 "figures/drug/Fig7_panel_TLR1_TIR_2D.pdf", is_pdf = TRUE)
# (elided: SLC25A46 3D/2D render panels and remaining fixed-size Fig7 panels)

# --- Composite layout (data panels) ---
theme_fig7 <- theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.3),
        axis.ticks.length = unit(1.5, "pt"),
        plot.title = element_blank(),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

NPG_TLR1 <- c("LRR" = "#C65D3B", "TIR" = "#2F6B5F")
NPG_SLC  <- c("Full" = "#2F6B5F", "MCF" = "#5BAF9A")
NPG_MD   <- c("TLR1 TIR" = "#2F6B5F", "SLC25A46" = "#C65D3B")

# Panel a: screening cascade
panel_a <- function() {
  d <- data.frame(
    stage = factor(c("ZINC\nPin1", "Unique\nSMILES", "Lipinski\n+PAINS",
                     "GraphBAN\nscored", "ADMET\npassed", "Vina\ndocked",
                     "MD\nvalidated")),
    count = c(498910, 249455, 237326, 237326, 164, 20, 2))
  ggplot(d, aes(x = stage, y = count)) +
    geom_col(fill = "#2F6B5F", width = 0.55) +
    geom_text(aes(label = comma(count)), vjust = -0.6, size = 2.5) +
    scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.22))) +
    theme_fig7 + labs(x = "", y = "Compounds (n)")
}

# Panels b/c: dual-domain docking comparison (top 10 per domain);
# single shared implementation, called once per target
panel_docking_compare <- function(target, domains, order_domain, hline, palette) {
  d <- rbindlist(lapply(names(domains), function(dm) {
    f <- file.path(DRUG_RESULTS, "docking", target, domains[[dm]], "docking_summary.csv")
    if (!file.exists(f)) return(NULL)
    dt <- fread(f); dt[, domain := dm]; dt
  }))
  d[, rank := frank(best_affinity_kcal_mol), by = domain]
  d <- d[rank <= 10]
  d[, lbl := factor(sub("rank_", "", ligand_id),
                    levels = d[domain == order_domain][order(-best_affinity_kcal_mol), sub("rank_", "", ligand_id)])]
  ggplot(d, aes(x = lbl, y = best_affinity_kcal_mol, fill = domain)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_hline(yintercept = hline, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_flip() + theme_fig7 +
    theme(legend.position = c(0.78, 0.2), legend.direction = "vertical",
          legend.key.size = unit(1.5, "mm"), legend.title = element_blank(),
          axis.text.y = element_text(size = 6)) +
    labs(x = "", y = "Vina affinity (kcal/mol)")
}
panel_b <- function()
  panel_docking_compare("TLR1", c("LRR domain" = "vina_top20_LRR",
                                  "TIR domain" = "vina_top20_TIR"),
                        "TIR", -6.0, NPG_TLR1)
panel_c <- function()
  panel_docking_compare("SLC25A46", c("Full domain" = "vina_top20_full",
                                      "MCF domain" = "vina_top20_MCF"),
                        "Full domain", -7.0, NPG_SLC)

# Panel d: MD backbone RMSD time series (2 lead systems)
panel_d <- function() {
  rmsd_all <- read_rmsd(MD_LEAD_MAP)
  ggplot(rmsd_all, aes(x = time_ns, y = rmsd_nm, color = system)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = NPG_MD) +
    scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_fig7 +
    theme(legend.position = c(0.8, 0.15), legend.direction = "vertical",
          legend.key.size = unit(1.5, "mm"), legend.title = element_blank()) +
    labs(x = "Time (ns)", y = "Backbone RMSD (nm)")
}

# Panel e: MD stability metrics (dot-plot comparison of the 2 lead systems)
panel_e <- function() {
  summaries <- lapply(MD_LEAD_MAP, function(sys) {
    f <- file.path(DRUG_RESULTS, "md", sys, "analysis_basic", "analysis_summary.json")
    if (file.exists(f)) fromJSON(f) else NULL
  })
  gm <- function(s, path) { val <- s[[path]]$mean; if (is.null(val)) NA_real_ else as.numeric(val) }
  metrics <- list(
    list(name = "RMSD (nm)",          path = "backbone_rmsd_nm",              scale = 1),
    list(name = "RMSF (nm)",          path = "ca_rmsf_nm",                    scale = 1),
    list(name = "Rg (nm)",            path = "protein_rg_nm",                 scale = 1),
    list(name = "Mindist (nm)",       path = "protein_ligand_mindist_nm",     scale = 1),
    list(name = "Contacts (per 100)", path = "protein_ligand_contacts_0p45nm", scale = 100))
  d <- melt(rbindlist(lapply(metrics, function(m) data.frame(
        metric = m$name,
        v1 = gm(summaries[[1]], m$path) / m$scale,
        v2 = gm(summaries[[2]], m$path) / m$scale))),
      id.vars = "metric", variable.name = "system", value.name = "value")
  d[, system := ifelse(system == "v1", names(MD_LEAD_MAP)[1], names(MD_LEAD_MAP)[2])]
  d[, metric := factor(metric, levels = sapply(metrics, `[[`, "name"))]
  ggplot(d, aes(x = value, y = metric, color = system)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = NPG_MD) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    theme_fig7 + labs(x = "", y = "")
}

# --- Composite assembly ---
p_a <- panel_a() + labs(tag = "a")
p_b <- panel_b() + labs(tag = "b")
p_c <- panel_c() + labs(tag = "c")
p_d <- panel_d() + labs(tag = "d")
p_e <- panel_e() + labs(tag = "e")

design <- "
AAAA
BBCC
DDEE
"

fig7 <- wrap_plots(A = p_a, B = p_b, C = p_c, D = p_d, E = p_e, design = design) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

ggsave(file.path(DRUG_FIGURES, "Fig7.pdf"), fig7,
       width = 180, height = 200, units = "mm", device = "pdf")

# -----------------------------------------------------------------------------
# 7. Supplementary tables (Table S13)
# -----------------------------------------------------------------------------
# Candidate rankings per target (ADMET + docking columns merged) and the MD
# summary table; written as tab-separated files under tables/drug/.

write_candidate_table <- function(target, preferred_domain, output_file) {
  admet_f <- file.path(DRUG_RESULTS, "admet", target, paste0(target, "_admet_passed.csv"))
  vina_f <- file.path(DRUG_RESULTS, "docking", target, preferred_domain, "docking_summary.csv")
  if (!file.exists(admet_f)) { warning("ADMET file not found: ", admet_f); return() }
  admet <- fread(admet_f)
  out_cols <- intersect(
    c("rank_after_admet", "SMILES", "graphban_mean_score",
      "biosnap_score", "kiba_score",
      "logP", "Molecular_Weight", "HBA", "HBD", "TPSA",
      "AMES", "Carcinogens_Lagunin", "ClinTox", "DILI", "hERG",
      "Bioavailability_Ma", "Half_Life_Obach", "PPBR_AZ"),
    names(admet))
  if (file.exists(vina_f)) {
    vina <- fread(vina_f)
    admet <- merge(admet, vina[, .(ligand_id, best_affinity_kcal_mol)],
                   by.x = "SMILES", by.y = "ligand_id", all.x = TRUE)
    out_cols <- c(out_cols, "best_affinity_kcal_mol")
  }
  out_cols <- intersect(out_cols, names(admet))
  fwrite(admet[, ..out_cols][order(rank_after_admet)], output_file, sep = "\t")
  message("Wrote ", basename(output_file))
}

write_md_summary_table <- function(output_file) {
  md_summary <- data.frame(
    System = c("TLR1 LRR", "TLR1 TIR", "SLC25A46 Full", "SLC25A46 MCF"),
    Vina_affinity = rep(NA_real_, 4),
    RMSD_mean_nm = rep(NA_real_, 4),
    RMSD_final_nm = rep(NA_real_, 4),
    Rg_final_nm = rep(NA_real_, 4),
    SASA_final_nm2 = rep(NA_real_, 4),
    PL_mindist_mean_nm = rep(NA_real_, 4),
    PL_contacts_mean = rep(NA_real_, 4),
    Hbonds_mean = rep(NA_real_, 4),
    MMGBSA_dG = rep(NA_real_, 4))
  sys_dirs <- c("TLR1_LRR_rank001", "TLR1_TIR_rank001",
                "SLC25A46_full_rank001", "SLC25A46_MCF_rank001")
  for (i in seq_along(sys_dirs)) {
    f <- file.path(DRUG_RESULTS, "md", sys_dirs[i], "analysis_basic", "analysis_summary.json")
    if (!file.exists(f)) next
    s <- fromJSON(f)
    md_summary$RMSD_mean_nm[i]   <- s$backbone_rmsd_nm$mean
    md_summary$RMSD_final_nm[i]  <- s$backbone_rmsd_nm$final
    md_summary$Rg_final_nm[i]    <- s$protein_rg_nm$final
    md_summary$SASA_final_nm2[i] <- s$protein_sasa_nm2$final
    md_summary$PL_mindist_mean_nm[i] <- s$protein_ligand_mindist_nm$mean
    md_summary$PL_contacts_mean[i]   <- s$protein_ligand_contacts_0p45nm$mean
    md_summary$Hbonds_mean[i] <- s$protein_ligand_hbond$mean_hbonds_per_frame
  }
  fwrite(md_summary, output_file, sep = "\t")
  message("MD summary table written: ", output_file)
}

generate_all_tables <- function() {
  write_candidate_table("TLR1", "vina_top20_LRR",
    file.path(DRUG_TABLES, "Table_drug_TLR1_candidates.tsv"))
  write_candidate_table("SLC25A46", "vina_top20_MCF",
    file.path(DRUG_TABLES, "Table_drug_SLC25A46_candidates.tsv"))
  write_md_summary_table(file.path(DRUG_TABLES, "Table_drug_MD_summary.tsv"))
}
