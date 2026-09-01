# =============================================================================
# 05_scrnaseq_tcr.R
# Results 2.5: Single-cell and paired TCR analyses linked prioritized genes to
# baseline response state and repertoire remodeling (paper Figure 6,
# Supplementary Figures S1-S10)
# Paper: Genetic prioritization of CD4+ T cell-expressed candidate genes for
#        allergic rhinitis (Biomedicines)
# Merged and sanitized from: sc_05-sc_06 (load/QC/integration of GSE200107,
#        shown in full; GSE273975 pipeline elided), sc_09-sc_11 (annotation,
#        MR landscape, pseudobulk),
#        sc_12-sc_21 (supplementary modules), sc_22-sc_33 (Ro/e, Fisher,
#        miloR, final figures), plot_density_v5.R
# =============================================================================
# Notes:
#  - Run from the repository root; each script assumes cwd == repo root.
#  - Raw scRNA-seq data are obtained from public repositories (GEO GSE200107
#    and GSE273975; see the paper's Data Availability statement); they are
#    not distributed here.
#  - Shared helpers live in scripts/R_functions_roe_analysis.R.
#  - Repetitive blocks are elided and marked "# (elided: ...)".
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries, paths and shared constants
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(stringr)
  library(scales)
  library(lme4)            # pseudobulk mixed models
  library(future)
  library(future.apply)
})

# Shared helpers (sourced, not copied): calc_gene_fisher, run_milor_da,
# plot_milor_beeswarm, load_target_genes, theme_nature, NPG_* colours.
source("scripts/R_functions_roe_analysis.R")

# Repo-relative paths (raw data: GEO; results/figures/tables: this repo)
RAW  <- "data"
RES  <- "results"
FIG  <- "figures"
SCRN <- file.path(RES, "scRNA", "GSE200107")

for (d in c(FIG, SCRN)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# -----------------------------------------------------------------------------
# 1. QC + integration, GSE200107 (representative dataset, shown in full)
#    GSE200107: CD4+ memory T cells, 15 libraries from 7 Japanese cedar
#    pollinosis patients (14 paired pre/post SLIT libraries + 1 extra post
#    library 46B-2); 4 good responders (25, 44, 46, 53) vs 3 poor (20, 48, 50).
# -----------------------------------------------------------------------------

# -- 1.1 Sample table (GEO GSM6008809-GSM6008823) ---------------------------
samples <- data.frame(
  gsm       = c("GSM6008809","GSM6008810","GSM6008811","GSM6008812",
                "GSM6008813","GSM6008814","GSM6008815","GSM6008816",
                "GSM6008817","GSM6008818","GSM6008819","GSM6008820",
                "GSM6008821","GSM6008822","GSM6008823"),
  sample_id = c("20A","20B","48A","48B","50A","50B","25A","25B","44A","44B","46A","46B","46B-2","53A","53B"),
  patient   = c("20","20","48","48","50","50","25","25","44","44","46","46","46","53","53"),
  condition = c("Pre","Post","Pre","Post","Pre","Post","Pre","Post","Pre","Post","Pre","Post","Post","Pre","Post"),
  stringsAsFactors = FALSE
)

# -- 1.2 Load 10X matrices and merge -----------------------------------------
seu_list <- list()
for (i in seq_len(nrow(samples))) {
  prefix <- paste0(samples$gsm[i], "_", samples$sample_id[i], "_")
  mat  <- readMM(file.path(RAW, "scRNA_GSE200107_RAW", paste0(prefix, "matrix.mtx.gz")))
  barcodes <- read.table(file.path(RAW, "scRNA_GSE200107_RAW", paste0(prefix, "barcodes.tsv.gz")),
                         header = FALSE, stringsAsFactors = FALSE)[, 1]
  features <- read.table(file.path(RAW, "scRNA_GSE200107_RAW", paste0(prefix, "features.tsv.gz")),
                         header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) <- paste0(samples$sample_id[i], "_", barcodes)
  rownames(mat) <- make.unique(features[, 2])

  obj <- CreateSeuratObject(counts = mat, project = samples$sample_id[i],
                            min.cells = 0, min.features = 0)
  obj$sample_id <- samples$sample_id[i]
  obj$patient   <- samples$patient[i]
  obj$condition <- samples$condition[i]
  seu_list[[i]] <- obj
}
merged <- merge(seu_list[[1]], seu_list[-1], project = "GSE200107")
merged <- JoinLayers(merged)

# -- 1.3 QC metrics and filters ----------------------------------------------
# QC thresholds (as used for the paper): genes >= 300, counts >= 500,
# MT% < 20, per-sample median absolute deviation 5x, gene min_cells = 10.
mt_genes <- grep("^MT-", rownames(merged), value = TRUE)
merged$pct_MT <- Matrix::colSums(merged[["RNA"]]$counts[mt_genes, , drop = FALSE]) /
                 Matrix::colSums(merged[["RNA"]]$counts) * 100
merged$n_genes_by_counts <- Matrix::colSums(merged[["RNA"]]$counts > 0)
merged$total_counts      <- Matrix::colSums(merged[["RNA"]]$counts)

keep <- merged$n_genes_by_counts >= 300 &
        merged$total_counts >= 500 &
        merged$pct_MT < 20

cells_keep <- c()
for (s in unique(merged$sample_id)) {
  idx <- which(merged$sample_id == s & keep)
  if (length(idx) == 0) next
  m_n <- merged$n_genes_by_counts[idx] >= median(merged$n_genes_by_counts[idx]) - 5 * mad(merged$n_genes_by_counts[idx], constant = 1) &
         merged$n_genes_by_counts[idx] <= median(merged$n_genes_by_counts[idx]) + 5 * mad(merged$n_genes_by_counts[idx], constant = 1)
  m_t <- merged$total_counts[idx] >= median(merged$total_counts[idx]) - 5 * mad(merged$total_counts[idx], constant = 1) &
         merged$total_counts[idx] <= median(merged$total_counts[idx]) + 5 * mad(merged$total_counts[idx], constant = 1)
  cells_keep <- c(cells_keep, idx[m_n & m_t])
}
merged <- merged[, cells_keep]
merged <- merged[Matrix::rowSums(merged[["RNA"]]$counts > 0) >= 10, ]

saveRDS(merged, file.path(SCRN, "01_seurat_qc.rds"))

# -- 1.4 TCR metadata integration (sc_06 Part A) -----------------------------
# High-confidence productive TRA/TRB chains per cell, collapsed to one row
# per barcode; clonotype_id from Cell Ranger raw_clonotype_id.
tcr_files <- list.files(file.path(RAW, "scRNA_GSE200107_RAW"),
                        pattern = "_all_contig_annotations.csv.gz$")
tcr_smp <- sapply(strsplit(tcr_files, "_"), function(x) x[2])
names(tcr_files) <- tcr_smp

for (col in c("TRA_v", "TRA_j", "TRA_cdr3", "TRB_v", "TRB_j", "TRB_cdr3")) merged[[col]] <- NA_character_
merged$clonotype_id <- NA_character_
merged$has_TCR <- FALSE

for (smp in names(tcr_files)) {
  tcr_raw <- read.csv(gzfile(file.path(RAW, "scRNA_GSE200107_RAW", tcr_files[smp])),
                      stringsAsFactors = FALSE)
  tcr_sub <- tcr_raw[tcr_raw$high_confidence == "true" & tcr_raw$productive == "true", ]
  if (nrow(tcr_sub) == 0) next

  collapse_chains <- function(df) {
    if (nrow(df) == 0) return(data.frame(barcode = character(), v = character(),
                                         j = character(), cdr3 = character(),
                                         stringsAsFactors = FALSE))
    agg <- aggregate(df[, c("v_gene", "j_gene", "cdr3")],
                     by = list(barcode = df$barcode),
                     FUN = function(x) paste(unique(x), collapse = "; "))
    colnames(agg)[2:4] <- c("v", "j", "cdr3")
    agg
  }
  tcr_tab <- merge(collapse_chains(tcr_sub[tcr_sub$chain == "TRA", ]),
                   collapse_chains(tcr_sub[tcr_sub$chain == "TRB", ]),
                   by = "barcode", all = TRUE, suffixes = c("_TRA", "_TRB"))
  bc_ct <- aggregate(raw_clonotype_id ~ barcode,
                     data = tcr_sub[tcr_sub$raw_clonotype_id != "", , drop = FALSE],
                     FUN = function(x) x[1])
  tcr_tab <- merge(tcr_tab, bc_ct, by = "barcode", all.x = TRUE)
  colnames(tcr_tab)[ncol(tcr_tab)] <- "clonotype_id"

  tcr_tab$cell_name <- paste0(smp, "_", tcr_tab$barcode)
  tcr_tab <- tcr_tab[tcr_tab$cell_name %in% colnames(merged), ]
  if (nrow(tcr_tab) == 0) next

  cell_idx <- match(tcr_tab$cell_name, colnames(merged))
  merged$has_TCR[cell_idx] <- TRUE
  for (col in c("v_TRA", "j_TRA", "cdr3_TRA", "v_TRB", "j_TRB", "cdr3_TRB", "clonotype_id")) {
    target <- c(v_TRA = "TRA_v", j_TRA = "TRA_j", cdr3_TRA = "TRA_cdr3",
                v_TRB = "TRB_v", j_TRB = "TRB_j", cdr3_TRB = "TRB_cdr3",
                clonotype_id = "clonotype_id")[[col]]
    vals <- tcr_tab[[col]]
    vals[is.na(vals) | vals == ""] <- NA_character_
    merged[[target]][cell_idx] <- vals
  }
}

# -- 1.5 Harmony re-integration + clustering (sc_09) --------------------------
# Initial LogNormalize pass was superseded by patient-aware Harmony
# re-integration; clustering evaluated at resolutions 0.2-1.2.
obj <- merged
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample_id)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))

saveRDS(obj, file.path(SCRN, "03_seurat_reint_harmony.rds"))

# -- 1.6 PBMC reference dataset (GSE273975): identical pipeline ----------------
# (elided: GSE273975 [sc_01/sc_02] runs the identical load/QC/integration
#  pipeline used for its descriptive PBMC reference role in the paper;
#  QC thresholds: genes>200/MT<20%/HB<5%.)

# -----------------------------------------------------------------------------
# 2. Cell-type annotation (multi-method: SingleR + reference markers + manual)
# -----------------------------------------------------------------------------

# -- 2.1 GSE200107: marker-based manual annotation of 20 CD4+ T subtypes ------
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
Idents(obj) <- "seurat_clusters"

markers_all <- FindAllMarkers(obj, assay = "RNA", slot = "data",
                              only.pos = TRUE, min.pct = 0,
                              logfc.threshold = 0, test.use = "wilcox")
markers_strict <- markers_all %>%
  filter(p_val_adj < 0.01, avg_log2FC > 1, pct.1 > 0.5, pct.2 < 0.5)

write.csv(markers_all, file.path(SCRN, "04_01_FindAllMarkers_all.csv"), row.names = FALSE)
write.csv(markers_strict, file.path(SCRN, "04_02_FindAllMarkers_strict.csv"), row.names = FALSE)

# 20 annotated subtypes (clusters 17 and 21 removed as contaminants)
cluster_anno <- c(
  "0"  = "C01_Tcm_CCR7",           # central memory T (CCR7)
  "1"  = "C02_Th2_Early_CCR4",     # Th2-prone memory T (PASK+, CCR4+)
  "2"  = "C03_Tem_Th1-like",       # Th1-like effector memory (CCL5+, CXCR3+)
  "3"  = "C04_Tem_Activated",      # activated effector memory (CAPG+, LGALS1+)
  "4"  = "C05_Th17_CCR6",          # Th17 (CCR6+, TNFSF13B+)
  "5"  = "C06_T_ISG_STAT1",        # interferon-responsive T (IFN signature)
  "6"  = "C07_Tcyto_GZMK",         # GZMK+ cytotoxic CD4+ T
  "7"  = "C08_T_Activated_CISH",   # early-activated memory T (CISH+)
  "8"  = "C09_Th2_Effector_GATA3", # pathogenic effector Th2 (GATA3+)
  "9"  = "C10_T_Activated_CD7",    # CD7-high activated memory T
  "10" = "C11_T_Metabolic_Low",    # low-metabolic memory T (low ribosomal)
  "11" = "C12_T_Apoptotic_BAX",    # apoptosis-prone T (BAX+)
  "12" = "C13_Treg_Activated",     # activated effector Treg (HLA-DR+, FOXP3+)
  "13" = "C14_T_Quiescent",        # quiescent long-lived memory T (MTRNR2L12+)
  "14" = "C15_CTL_GNLY",           # terminally differentiated CTL-like CD4+ T
  "15" = "C16_T_Metabolic_High",   # high-metabolic T (high ribosomal)
  "16" = "C17_Trm-like_GZMK",      # tissue-resident-like memory T (GZMK)
  "18" = "C18_Treg_IKZF2",         # thymic/steady-state Treg (Helios/IKZF2+)
  "19" = "C19_T_Stressed_HIF1A",   # hypoxic-stress / transcriptionally active T
  "20" = "C20_T_Cycling"           # proliferating T (MKI67+, RRM2+)
)

obj <- subset(obj, idents = names(cluster_anno))
obj <- RenameIdents(obj, cluster_anno)
obj$cell_type <- Idents(obj)

# Patient-level response assignment (4 good vs 3 poor responders)
obj$response <- ifelse(obj$patient %in% c("20", "48", "50"),
                       "Poor Response", "Good Response")
obj$response <- factor(obj$response, levels = c("Good Response", "Poor Response"))

saveRDS(obj, file.path(SCRN, "05_manual_annotation.rds"))

# -- 2.2 GSE273975 annotation -------------------------------------------------
# (elided: identical workflow [sc_09_annotation_GSE273975.R] yielding 24 PBMC
#  types across 16,412 PBMCs / 4 libraries / 2 patients; descriptive
#  reference dataset, not a replication cohort.)

# -----------------------------------------------------------------------------
# 3. MR gene expression landscape in scRNA (GSE200107)
# -----------------------------------------------------------------------------
sce <- readRDS(file.path(SCRN, "05_manual_annotation.rds"))
DefaultAssay(sce) <- "RNA"
sce <- JoinLayers(sce, assay = "RNA")

# Tier label per gene (helper for dotplot ordering)
tier_group_fill <- function(mr_dt, tg) {
  lab <- rep("Other", nrow(mr_dt))
  for (nm in names(tg)[1:6]) {
    lab[mr_dt$gene_symbol %in% tg[[nm]]] <- nm
  }
  lab
}

mr_full <- fread(file.path(RES, "05_MR_ivw_fdr.tsv"))
sym_dt  <- fread(file.path(RES, "gene_symbol_full.tsv"))
setnames(sym_dt, c("gene", "gene_symbol"))
sym_lup <- setNames(sym_dt$gene_symbol, sym_dt$gene)

mr_best <- mr_full[fdr < 0.05][order(fdr, pval)][!duplicated(gene)]
mr_best[, gene_symbol := coalesce(sym_lup[gene], gene)]
mr_best[, direction := ifelse(beta > 0, "Risk", "Protect")]
mr_genes <- mr_best[gene_symbol %in% rownames(sce)]

tier_genes <- list(
  Tier1_Risk    = c("TSLP", "ORMDL3"),
  Tier1_Protect = c("SLC25A46", "TLR1"),
  Tier2_Risk    = c("IL18R1", "IL18RAP"),
  Tier2_Protect = c("IKZF3", "CD247"),
  Tier3_Risk    = c("D2HGDH", "STAT6", "LMBRD2"),
  Tier3_Protect = c("PDCD1", "IL1R2")
)
tier_genes$All_Risk    <- mr_genes[direction == "Risk", gene_symbol]
tier_genes$All_Protect <- mr_genes[direction == "Protect", gene_symbol]
tier_genes <- lapply(tier_genes, function(g) intersect(g, rownames(sce)))
tier_genes <- tier_genes[lengths(tier_genes) > 0]

# Module scoring (AddModuleScore appends "1" to the supplied name)
for (nm in names(tier_genes)) {
  col_name <- paste0("score_", nm)
  sce <- AddModuleScore(sce, features = list(tier_genes[[nm]]),
                        name = col_name, ctrl = min(100, nrow(sce)))
  sce@meta.data[[col_name]] <- sce@meta.data[[paste0(col_name, "1")]]
  sce@meta.data[[paste0(col_name, "1")]] <- NULL
}

# MR gene DotPlot (Supp Fig S5)
dot_genes_tier <- mr_genes[tier_group_fill(mr_genes, tier_genes) != "Other", gene_symbol]
p_dot <- DotPlot(sce, features = dot_genes_tier,
                 group.by = "cell_type", split.by = "condition",
                 cols = c("lightgrey", NPG_RED)) +
  RotatedAxis() +
  labs(title = "MR Tier 1-3 gene expression across CD4+ T subclusters",
       subtitle = "Split by condition (Pre / Post SLIT)") +
  theme_nature(9) +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom")
ggsave(file.path(FIG, "scRNA_MR_dotplot_Tier_genes.pdf"),
       p_dot, width = max(10, length(dot_genes_tier) * 0.35), height = 6,
       limitsize = FALSE)

# (elided: per-module trend plots, module-score heatmap and Layer-2 candidate
#  screening from sc_10; outputs: 06_module_scores.csv, 07_layer2_candidates.tsv.)
# (elided: GSE273975 MR landscape [sc_15]: 115/129 MR genes detected in PBMCs,
#  TLR1 monocyte-enriched (Ro/e = 2.52), CD247 highest positive rate (58%).)
# (elided: novel-gene characterization [sc_21] for LMBRD2, MARS2, KLHL5,
#  RTP5, SARNP; exploratory, without statistics.)

# -----------------------------------------------------------------------------
# 4. Positive-cell enrichment, baseline Poor-Pre vs Good-Pre (paper Fig 6c)
#    Ro/e (ratio of observed to expected) + Fisher exact test per cell type,
#    BH FDR. Cell-level result: 561 gene-cell-type pairs pass FDR < 0.05.
# -----------------------------------------------------------------------------
# (elided: Pre-vs-Post overall Ro/e + Fisher [sc_22] and Good-vs-Poor overall
#  Ro/e [sc_24 analyses A-C]; superseded by the paired analysis in Section 5.
#  Identical parameters: positive = expression > 0, Fisher 2x2 within cell
#  type, min_pos_total = 5, min_pos_per_group = 3, BH FDR < 0.05.)

gene_dt <- load_target_genes(".")
valid_genes <- gene_dt[gene_symbol %in% rownames(sce)]
gene_list <- valid_genes$gene_symbol
tier123_genes <- valid_genes[tier %in% c("Tier1", "Tier2", "Tier3")]$gene_symbol

# Pre-extract expression (avoids shipping the Seurat object to workers)
expr_all <- FetchData(sce, vars = c(gene_list, "cell_type", "condition",
                                    "patient", "response"))
expr_all$patient <- as.character(expr_all$patient)

# Baseline contrast: Poor-Pre vs Good-Pre (only Pre cells). As in sc_26,
# cond1 = Poor Pre cells (first response group encountered in cell order);
# log2_roe_ratio > 0 therefore means Poor-enriched (red in the heatmap).
baseline_df <- expr_all[expr_all$condition == "Pre", ]

n_cores <- max(1, parallel::detectCores() - 2)
plan("multisession", workers = n_cores)

baseline_fisher <- future_lapply(seq_along(gene_list), function(i) {
  calc_gene_fisher(expr_df = baseline_df, gene = gene_list[i],
                   celltype_col = "cell_type", condition_col = "response",
                   min_pos_total = 5, min_pos_per_group = 3)
}, future.seed = 42L)

bl_dt <- rbindlist(baseline_fisher, fill = TRUE)
bl_dt[, fdr := p.adjust(fisher_pval, method = "BH")]
bl_dt[, significant := fdr < 0.05]
bl_dt[, enrich_group := fifelse(log2_roe_ratio > 0, "Poor Response",
                         fifelse(log2_roe_ratio < 0, "Good Response", "NS"))]
bl_dt <- merge(bl_dt, valid_genes[, .(gene_symbol, tier, direction)],
               by.x = "gene", by.y = "gene_symbol", all.x = TRUE)
setnames(bl_dt,
         c("n_cond1", "n_cond2", "n_cond1_pos", "n_cond1_neg",
           "n_cond2_pos", "n_cond2_neg", "pct_cond1_pos", "pct_cond2_pos",
           "roe_cond1", "roe_cond2"),
         c("n_Poor_Pre", "n_Good_Pre", "n_Poor_Pre_pos", "n_Poor_Pre_neg",
           "n_Good_Pre_pos", "n_Good_Pre_neg", "pct_Poor_Pre", "pct_Good_Pre",
           "roe_Poor_Pre", "roe_Good_Pre"))

fwrite(bl_dt, file.path(SCRN, "10_baseline_poor_vs_good.tsv"), sep = "\t")
fwrite(bl_dt[fdr < 0.05][order(fdr)], file.path(SCRN, "10_baseline_poor_vs_good_fdr.tsv"), sep = "\t")

# Result (paper Results 2.5): 561 significant gene-cell-type pairs at cell
# level (FDR < 0.05, 88 unique genes). Protective genes SLC25A46, TLR1,
# CD247, IKZF3 are consistently LOWER in Poor responders; RPS26 is
# systematically higher in Poor; ORMDL3 (risk gene) is paradoxically higher
# in Good responders. CD247 is higher in Good across 14/20 subtypes (e.g.
# C05_Th17_CCR6 57.3% vs 40.1%; C01_Tcm_CCR7 54.2% vs 43.1%); TLR1 is
# largely confined to activated compartments (C10_T_Activated_CD7 11.8% vs 5.7%).

# -----------------------------------------------------------------------------
# 5. Paired Pre-to-Post delta analysis (paper Fig 6d)
#    Per-patient positive-cell proportions per gene-cell-type pair;
#    delta = Post - Pre; Wilcoxon test Poor delta vs Good delta (BH FDR).
# -----------------------------------------------------------------------------
prop_list <- future_lapply(seq_along(gene_list), function(i) {
  sub <- expr_all[, c(gene_list[i], "cell_type", "condition", "patient", "response")]
  names(sub)[1] <- "expr"
  sub$positive <- sub$expr > 0
  sub <- as.data.table(sub)
  agg <- sub[, .(n_cells = .N, n_pos = sum(positive), pct_pos = mean(positive) * 100),
             by = .(patient, cell_type, condition)]
  agg$gene <- gene_list[i]
  merge(agg, unique(sub[, .(patient, response)]), by = "patient")
}, future.seed = 42L)
plan("sequential")

prop_dt <- rbindlist(prop_list, fill = TRUE)
prop_wide <- dcast(prop_dt, gene + patient + cell_type + response ~ condition,
                   value.var = "pct_pos", fill = NA)
prop_wide[, delta := Post - Pre]
prop_wide <- prop_wide[!is.na(delta)]

# 15,827 paired observations (gene x patient x cell type)
fwrite(prop_wide, file.path(SCRN, "10_paired_delta.tsv"), sep = "\t")

delta_test <- prop_wide[, {
  poor_delta <- delta[response == "Poor Response"]
  good_delta <- delta[response == "Good Response"]
  if (length(poor_delta) >= 2 && length(good_delta) >= 2) {
    wt <- tryCatch(wilcox.test(poor_delta, good_delta, exact = FALSE), error = function(e) NULL)
    .(n_poor = length(poor_delta), n_good = length(good_delta),
      mean_delta_poor = mean(poor_delta, na.rm = TRUE),
      mean_delta_good = mean(good_delta, na.rm = TRUE),
      wilcox_pval = if (!is.null(wt)) wt$p.value else NA_real_)
  } else {
    .(n_poor = length(poor_delta), n_good = length(good_delta),
      mean_delta_poor = NA_real_, mean_delta_good = NA_real_, wilcox_pval = NA_real_)
  }
}, by = .(gene, cell_type)]

delta_test <- merge(delta_test, valid_genes[, .(gene_symbol, tier, direction)],
                    by.x = "gene", by.y = "gene_symbol", all.x = TRUE)
delta_test[, fdr_wilcox := p.adjust(wilcox_pval, method = "BH")]
delta_test[, significant := fdr_wilcox < 0.05]
fwrite(delta_test, file.path(SCRN, "10_delta_wilcoxon.tsv"), sep = "\t")

# Result (paper Results 2.5): 0 gene-cell-type pairs survive FDR at the
# patient level (8 nominal p < 0.05, ~113 expected by chance); treatment
# slopes are parallel between groups - baselines differ, changes do not.

# -- 5.1 Paired-delta figure: top 8 gene-cell-type pairs (Fig 6d) -------------
RED  <- "#C62828"
BLUE <- "#1565C0"

candidates <- bl_dt[significant == TRUE & cell_type != ""]
candidates[, `:=`(abs_lor = abs(log2_roe_ratio),
                  tier_weight = fifelse(tier == "Tier1", 3,
                                 fifelse(tier == "Tier2", 2,
                                 fifelse(tier == "Tier3", 1, 0))))]
setorder(candidates, -tier_weight, -abs_lor)
picked <- rbind(candidates[tier == "Tier1"][1:3],
                candidates[tier == "Tier2"][1:2],
                candidates[gene == "RPS26"][1],
                candidates[gene == "MYL6" & cell_type == "C14_T_Quiescent"],
                candidates[gene == "NFKB1"][1])
picked <- unique(picked, by = c("gene", "cell_type"))[1:8]

delta_pairs <- prop_wide[gene %in% picked$gene & cell_type %in% picked$cell_type]
delta_long <- melt(delta_pairs, id.vars = c("gene", "patient", "cell_type", "response"),
                   measure.vars = c("Pre", "Post"),
                   variable.name = "timepoint", value.name = "pct_pos")
delta_long[, response := factor(response, levels = c("Poor Response", "Good Response"))]
delta_long <- merge(delta_long, unique(picked[, .(gene, cell_type, tier, direction)]), by = c("gene", "cell_type"))
panel_order <- unique(delta_long[, .(panel = paste0(gene, "\n", cell_type), tier, direction)])
panel_order <- panel_order[order(tier, direction)]
delta_long[, panel := factor(paste0(gene, "\n", cell_type), levels = panel_order$panel)]

p_delta <- ggplot(delta_long, aes(x = timepoint, y = pct_pos, group = patient)) +
  geom_line(aes(colour = response), linewidth = 0.6, alpha = 0.7) +
  geom_point(aes(colour = response), size = 2.5, alpha = 0.9) +
  stat_summary(aes(group = response, colour = response), fun = mean,
               geom = "line", linewidth = 1.8, alpha = 0.9) +
  stat_summary(aes(group = response, colour = response), fun = mean,
               geom = "point", size = 4, shape = 18) +
  scale_colour_manual(values = c("Poor Response" = RED, "Good Response" = BLUE),
                      name = NULL) +
  scale_x_discrete(expand = c(0.15, 0.15)) +
  facet_wrap(~ panel, ncol = 4, scales = "free_y") +
  labs(title = "Per-patient Pre-to-Post change in % positive cells",
       subtitle = paste0("Red = Poor Response (n=3), Blue = Good Response (n=4); ",
                         "thin lines = individual patients, thick diamonds = group mean; ",
                         "baselines differ, treatment slopes are similar."),
       x = NULL, y = "% positive cells") +
  theme_nature(9) + theme(legend.position = "top",
                          strip.text = element_text(size = 7, face = "bold"))
ggsave(file.path(FIG, "Fig6d_paired_delta.pdf"), p_delta,
       width = 320, height = 180, units = "mm", limitsize = FALSE)

# -----------------------------------------------------------------------------
# 6. Patient-level pseudobulk validation (mixed model)
#    Per patient x condition pseudobulk means;
#    lmer(mean_expr ~ response * condition + (1 | patient)); Wald p; BH FDR.
# -----------------------------------------------------------------------------
candidates_pb <- fread(file.path(SCRN, "07_layer2_candidates.tsv"))
if (nrow(candidates_pb) == 0) {
  # Fallback: all Tier 1+2 genes x all cell types
  tier12 <- c("TSLP", "ORMDL3", "SLC25A46", "TLR1", "IL18R1", "IL18RAP", "IKZF3", "CD247")
  candidates_pb <- as.data.table(expand.grid(gene = tier12,
                                             cell_type = levels(sce$cell_type),
                                             stringsAsFactors = FALSE))
}

pb_results <- list()
for (i in seq_len(nrow(candidates_pb))) {
  gene_i <- candidates_pb$gene[i]
  ct_i   <- candidates_pb$cell_type[i]
  if (!gene_i %in% rownames(sce)) next

  cells_ct <- WhichCells(sce, idents = ct_i)
  if (length(cells_ct) < 20) next
  sce_ct <- subset(sce, cells = cells_ct)

  expr_vec <- FetchData(sce_ct, vars = c(gene_i, "patient", "condition", "response"))
  colnames(expr_vec)[1] <- "expr"
  pb <- expr_vec %>%
    group_by(patient, condition, response) %>%
    summarise(mean_expr = mean(expr, na.rm = TRUE),
              n_cells = n(), .groups = "drop")
  if (any(pb$n_cells < 10) || nrow(pb) < 10) next

  pb$condition <- factor(pb$condition, levels = c("Pre", "Post"))
  pb$response  <- factor(pb$response,  levels = c("Good Response", "Poor Response"))
  pb$patient   <- factor(pb$patient)

  fit <- tryCatch(lmer(mean_expr ~ response * condition + (1 | patient), data = pb),
                  error = function(e) NULL)
  if (is.null(fit)) next

  coef_tab <- summary(fit)$coefficients
  int_row <- which(grepl("response.*condition|condition.*response", rownames(coef_tab)))
  if (length(int_row) == 0) next

  pb_results[[length(pb_results) + 1]] <- data.frame(
    gene = gene_i, cell_type = ct_i,
    n_patients = length(unique(pb$patient)), n_obs = nrow(pb),
    beta_interaction = coef_tab[int_row, "Estimate"],
    se_interaction   = coef_tab[int_row, "Std. Error"],
    p_interaction    = 2 * pnorm(-abs(coef_tab[int_row, "t value"])),
    delta_Good = mean(pb$mean_expr[pb$response == "Good Response" & pb$condition == "Post"]) -
                 mean(pb$mean_expr[pb$response == "Good Response" & pb$condition == "Pre"]),
    delta_Poor = mean(pb$mean_expr[pb$response == "Poor Response" & pb$condition == "Post"]) -
                 mean(pb$mean_expr[pb$response == "Poor Response" & pb$condition == "Pre"]),
    stringsAsFactors = FALSE)
}

res_dt <- rbindlist(pb_results)
res_dt[, p_fdr := p.adjust(p_interaction, method = "BH")]
res_dt <- res_dt[order(p_interaction)]
fwrite(res_dt, file.path(SCRN, "08_pseudobulk_results.tsv"), sep = "\t")

# Result (paper Results 2.5): 4 good vs 3 poor responders; minimum FDR = 0.209,
# two pairs reach nominal p < 0.05; no pair survives patient-level FDR.

# -----------------------------------------------------------------------------
# 7. miloR neighborhood differential abundance, Poor vs Good (paper Fig 6f)
#    k = 30, d = 20, prop = 0.2, SpatialFDR < 0.05 (helper: run_milor_da).
# -----------------------------------------------------------------------------
milo_res <- run_milor_da(sce, condition_col = "response", sample_col = "sample_id",
                         reduced_dim = "harmony", k = 30, d = 20, prop = 0.2)
da_results <- milo_res$da_results

# Result (paper Results 2.5): 5,384/17,789 neighborhoods DA at baseline
# (SpatialFDR < 0.05, 30.3%), strongest in C01_Tcm_CCR7 (2,050).
fwrite(as.data.table(da_results), file.path(SCRN, "11_milo_response_da.tsv"), sep = "\t")

# (elided: miloR Pre-vs-Post run with condition as contrast; 0/17,788
#  neighborhoods significant - treatment alone does not reshape
#  neighborhoods. Neighborhood graph plot also elided.)

p_milo <- plot_milor_beeswarm(da_results, group.by = "cell_type", alpha = 0.8) +
  labs(title = "miloR: Poor vs Good Response",
       subtitle = sprintf("%d significant neighborhoods (SpatialFDR < 0.05)", sum(da_results$SpatialFDR < 0.05, na.rm = TRUE)))
ggsave(file.path(FIG, "Fig6f_miloR_beeswarm.pdf"), p_milo, width = 12, height = 7)

# -----------------------------------------------------------------------------
# 8. Paired TCR clonotype dynamics (Supplementary Figures S1-S4)
#    Compact version of sc_13; 6 patients with matched pre/post TCR data.
#    Patient 53 is excluded from clonotype tracking (its post-treatment
#    library has no usable high-confidence TCR chains; pre library: 402 cells).
# -----------------------------------------------------------------------------
tcr_obj <- readRDS(file.path(SCRN, "05_manual_annotation.rds"))

response_map <- c("25" = "Good", "44" = "Good", "46" = "Good", "53" = "Good",
                  "20" = "Poor", "48" = "Poor", "50" = "Poor")
resp_vec <- unname(response_map[as.character(tcr_obj$patient)])
resp_vec[is.na(resp_vec)] <- "Unknown"
tcr_obj[["response_group"]] <- resp_vec

tcr_meta <- as.data.table(tcr_obj[[]], keep.rownames = "cell_barcode")
tcr_cells <- tcr_meta[has_TCR == TRUE]

# -- 8.1 Clonal diversity metrics (Shannon, inverse Simpson, Pielou, Gini) ----
compute_diversity <- function(dt) {
  clone_tab <- dt[, .N, by = .(clonotype_id)]
  total <- sum(clone_tab$N)
  freq <- clone_tab$N / total
  shannon <- -sum(freq * log(freq))
  inv_simpson <- 1 / sum(freq^2)
  n_clones <- nrow(clone_tab)
  sorted_n <- sort(clone_tab$N)
  gini <- (2 * sum(seq_len(n_clones) * sorted_n) / (n_clones * sum(sorted_n))) -
          (n_clones + 1) / n_clones
  data.table(n_cells = total, n_clones = n_clones,
             shannon = shannon, inv_simpson = inv_simpson,
             pielou = ifelse(n_clones > 1, shannon / log(n_clones), 0),
             gini = gini,
             expansion_pct = sum(clone_tab[N >= 2]$N) / total * 100)
}
diversity_metrics <- tcr_cells[, compute_diversity(.SD),
                               by = .(patient, condition, response_group)]
fwrite(diversity_metrics, file.path(SCRN, "TCR_diversity_metrics.tsv"), sep = "\t")
# (elided: diversity bar charts and per-metric faceted figures from sc_13;
#  no significant clonal diversity difference between groups, n = 7.)

# -- 8.2 Clone size groups and MR gene expression by expansion ----------------
clone_size_tot <- tcr_cells[, .N, by = .(clonotype_id)]
clone_size_tot[, expansion_group := fcase(N >= 5, "Expanded (>=5)",
                                          N >= 2 & N <= 4, "Small (2-4)",
                                          N == 1, "Singleton")]
tcr_cells <- merge(tcr_cells, clone_size_tot[, .(clonotype_id, expansion_group)],
                   by = "clonotype_id", all.x = TRUE)

core_genes <- valid_genes[core == TRUE]$gene_symbol   # 9 core genes
expr_mat <- GetAssayData(tcr_obj, assay = "RNA", layer = "data")[core_genes, tcr_cells$cell_barcode, drop = FALSE]
expr_long <- melt(as.data.table(as.matrix(expr_mat), keep.rownames = "gene"),
                  id.vars = "gene", variable.name = "cell_barcode", value.name = "expression")
expr_long <- merge(expr_long,
                   tcr_cells[, .(cell_barcode, patient, response_group, clonotype_id, expansion_group, cell_type)],
                   by = "cell_barcode", all.x = TRUE)

core_by_group <- expr_long[, .(pct_detect = sum(expression > 0) / .N * 100),
                           by = .(gene, expansion_group)]
# CD247 detected more often in expanded clones (44.9%) than singletons
# (40.5%) - the main TCR-MR gene association in the paper.
fwrite(core_by_group, file.path(SCRN, "TCR_core_genes_by_expansion.tsv"), sep = "\t")
# (elided: MR-gene heatmaps by expansion group, UMAP overlays and violins
#  from sc_13 figures 3/4/7/10/11; Supplementary Figures S5-S6 panels.)

# -- 8.3 Clonotype persistence and cell-type state transitions Pre -> Post ----
track_clonotypes <- function(pt_id) {
  pt_data <- tcr_cells[patient == pt_id]
  pre  <- unique(pt_data[condition == "Pre"]$clonotype_id)
  post <- unique(pt_data[condition == "Post"]$clonotype_id)
  list(shared = intersect(pre, post), pre_only = setdiff(pre, post),
       post_only = setdiff(post, pre), n_pre = length(pre), n_post = length(post))
}

transition_data <- rbindlist(lapply(unique(tcr_cells[patient != "53"]$patient), function(pt) {
  tr <- track_clonotypes(pt)
  resp <- unique(tcr_cells[patient == pt]$response_group)
  get_state <- function(cond, clones) {
    tcr_cells[patient == pt & condition == cond & clonotype_id %in% clones,
              .(dominant = names(sort(table(cell_type), decreasing = TRUE))[1],
                n_cells = .N), by = .(clonotype_id)]
  }
  pre <- get_state("Pre", tr$shared)
  post <- get_state("Post", tr$shared)
  setnames(pre,  c("dominant", "n_cells"), c("pre_dominant", "pre_n"))
  setnames(post, c("dominant", "n_cells"), c("post_dominant", "post_n"))
  m <- merge(pre, post, by = "clonotype_id")
  m[, `:=`(patient = pt, response_group = resp,
           switched = pre_dominant != post_dominant,
           expanded = (post_n - pre_n) > 0.5 * pmax(pre_n, 1))]
  m
}))

new_clones_summary <- rbindlist(lapply(unique(tcr_cells[patient != "53"]$patient), function(pt) {
  tr <- track_clonotypes(pt)
  data.table(patient = pt, response_group = unique(tcr_cells[patient == pt]$response_group),
             n_new_clones = length(tr$post_only),
             pct_new_of_post = length(tr$post_only) / tr$n_post * 100)
}))
fwrite(new_clones_summary, file.path(SCRN, "TCR_new_clones.tsv"), sep = "\t")

# Result (paper Results 2.5): post-only clonotypes 23.3% (Good) vs 3.5%
# (Poor) of post-treatment clones, n = 3 per group, with patient 46 driving
# the Good-responder values; shared-clonotype expansion 7.5% vs 1.6%;
# dominant cell-type state switching 88.7% vs 90.1%. Remodeling is dominated
# by newly emerged clonotypes in Good responders, while Poor responders
# show mostly persistent, size-stable clones.
# (elided: clone-size distribution bars, alluvial persistence diagrams,
#  transition heatmaps and clone-size scatter plots from sc_13 figures
#  1/2/5/6/9/12; these feed Supplementary Figures S1-S4.)

# -----------------------------------------------------------------------------
# 9. Supplementary modules (heavily elided; each feeds one Supp figure)
# -----------------------------------------------------------------------------

# -- 9.1 CellChat ligand-receptor analysis, GSE200107 (sc_12) ----------------
# Representative block (Good responders); Poor responders and the merged
# comparison are identical. Findings: MIF and CD99 dominate both groups;
# MHC-II signaling higher in Good responders; 7 MR genes map to CellChatDB.
# GSE273975 CellChat [sc_16] is elided identically. Supp Fig S8.
#
#   obj_cc <- readRDS(file.path(SCRN, "05_manual_annotation.rds"))
#   Idents(obj_cc) <- "cell_type"
#   cells_cc <- colnames(obj_cc)[obj_cc$response == "Good Response"]
#   sampled <- unlist(lapply(unique(obj_cc$cell_type), function(ct) {
#     cc <- cells_cc[obj_cc$cell_type[cells_cc] == ct]
#     if (length(cc) > 500) sample(cc, 500) else cc
#   }))
#   obj_cc <- subset(obj_cc, cells = sampled)
#   cc <- createCellChat(object = obj_cc, group.by = "cell_type")
#   cc@DB <- CellChatDB.human
#   cc <- subsetData(cc); cc <- identifyOverExpressedGenes(cc)
#   cc <- identifyOverExpressedInteractions(cc)
#   cc <- computeCommunProb(cc); cc <- computeCommunProbPathway(cc)
#   cc <- aggregateNet(cc); saveRDS(cc, file.path(SCRN, "09_cellchat_good.rds"))

# -- 9.2 TF regulon activity, GSE200107 (sc_14) ------------------------------
# Findings: NFKB1 regulon higher in Good Response (delta = -0.05, p = 3e-97);
# STAT6 higher in Poor (delta = -0.009, p = 4e-28); lineage validation across
# GSE273975 confirms IKZF3 in B/T cells, PDCD1 in T/NK, TLR1 in monocytes.
# Supp Fig S9. Representative block:
#
#   tf_genes <- c("STAT6", "IKZF3", "ZBTB38", "NFKB1", "PDCD1", "IL1R2", "TLR1")
#   tf_obj <- readRDS(file.path(SCRN, "05_manual_annotation.rds"))
#   tf_obj$group <- paste(tf_obj$response, tf_obj$condition, sep = "_")
#   # decoupleR ULM + VIPER regulon inference (dorothea) -> per-TF dotplots

# -- 9.3 Mitochondrial MR module and PDCD1 exhaustion -------------------------
# (elided: sc_19 Mito-MR module [SLC25A46, D2HGDH, MARS2, RPS26] and sc_20
#  PDCD1 exhaustion signature [PDCD1 + exhaustion markers, Pre vs Post;
#  higher in exhausted-like T cells of Poor responders]. Both follow the
#  AddModuleScore -> per-cell-type aggregate -> plot pattern. Supplementary
#  figures only; single-gene figures were dropped as over-fitted.)

# -----------------------------------------------------------------------------
# 10. Final publication figures (paper Figure 6a-f)
#     Supplementary panel sources are mapped in Section 10.7.
# -----------------------------------------------------------------------------
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(colour = "black", size = base_size - 0.5),
          axis.line = element_line(linewidth = 0.3),
          legend.title = element_text(size = base_size - 0.5, face = "bold"),
          legend.text = element_text(size = base_size - 1),
          strip.text = element_text(size = base_size - 0.5, face = "bold"),
          strip.background = element_rect(fill = "grey93", colour = NA),
          panel.grid = element_blank(),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          plot.subtitle = element_text(size = base_size - 0.5, colour = "grey40"),
          legend.position = "right")
}
WHITE <- "#F8F8F8"

# Compact version of plot_density_v5.R: per-gene UMAP density overlays
# (Seurat V5; layer defaults to "data"; used for MR-gene overlay panels).
plot_density_v5 <- function(object, features, layer = "data", reduction = "umap",
                            dims = c(1, 2), aspect_ratio = 1, xlab = NULL,
                            ylab = NULL, pal = "viridis", size = 1,
                            raster = FALSE) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  emb <- as.data.frame(Embeddings(object, reduction)[, dims])
  colnames(emb) <- c("d1", "d2")
  dat <- FetchData(object, vars = features, layer = layer)
  plots <- lapply(features, function(f) {
    p <- ggplot(cbind(emb, feature = dat[[f]]), aes(x = d1, y = d2, colour = feature)) +
      geom_point(shape = 16, size = size) +
      scale_colour_viridis_c(option = pal, name = "Expression") +
      coord_fixed(ratio = aspect_ratio) +
      xlab(xlab %||% "UMAP 1") + ylab(ylab %||% "UMAP 2") +
      ggtitle(f) + theme_classic()
    if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
      ggrastr::rasterise(p, dpi = 300)
    } else p
  })
  if (length(plots) == 1) plots[[1]] else wrap_plots(plots)
}

# -- 10.1 Fig 6a: GSE200107 UMAP colored by the 20 annotated subtypes --------
g_6a <- DimPlot(sce, group.by = "cell_type", label = FALSE, pt.size = 0.1) +
  theme_pub(8) +
  labs(title = "GSE200107: CD4+ memory T cells (117,771 cells)",
       subtitle = "20 annotated subtypes, 7 JCP patients (15 libraries)") +
  theme(legend.text = element_text(size = 5.5), legend.key.size = unit(0.3, "cm"))
ggsave(file.path(FIG, "Fig6a_umap_GSE200107.pdf"), g_6a,
       width = 180, height = 140, units = "mm", limitsize = FALSE)

# -- 10.2 Fig 6b: GSE273975 UMAP (descriptive reference dataset) -------------
# (elided: load results/scRNA/GSE273975/05_manual_annotation.rds and plot
#  DimPlot(group.by = "cell_type") as in 10.1, 24 PBMC types, 16,412 cells.)
# g_6b <- DimPlot(sce_273975, group.by = "cell_type", ...)
# ggsave(file.path(FIG, "Fig6b_umap_GSE273975.pdf"), g_6b, ...)

# -- 10.3 Fig 6c: baseline Tier 1-3 heatmap (sc_29 FigA code) -----------------
bl <- fread(file.path(SCRN, "10_baseline_poor_vs_good.tsv"))
bl[, log2OR := pmax(pmin(log2_roe_ratio, 2.5), -2.5)]
bl[is.na(log2OR), log2OR := 0]
bl[, star := fifelse(fdr < 0.001, "***", fifelse(fdr < 0.01, "**", fifelse(fdr < 0.05, "*", "")))]

pa <- bl[gene %in% tier123_genes & cell_type != ""]
pa[, tier := factor(tier, c("Tier1", "Tier2", "Tier3"))]
ct_order <- setdiff(sort(unique(bl$cell_type)), c("", NA))
pa[, cell_type := factor(cell_type, rev(ct_order))]
go <- pa[, .(m = mean(abs(log2OR))), by = .(gene, tier, direction)][order(tier, direction, -m)]
pa[, gene := factor(gene, levels = go$gene)]
pa[, tier_dir := factor(paste0(tier, " / ", direction),
     levels = unique(paste0(go$tier, " / ", go$direction)))]

g_6c <- ggplot(pa, aes(x = gene, y = cell_type)) +
  geom_tile(aes(fill = log2OR), colour = "grey88", linewidth = 0.35) +
  geom_text(aes(label = star), size = 3.5, colour = "grey20") +
  scale_fill_gradient2(low = BLUE, mid = WHITE, high = RED, midpoint = 0,
                       limits = c(-2.5, 2.5), oob = squish,
                       name = "Poor enriched     Good enriched",
                       guide = guide_colourbar(title.position = "top",
                                               title.hjust = 0.5,
                                               barwidth = 0.5, barheight = 4)) +
  facet_grid(. ~ tier_dir, scales = "free_x", space = "free_x") +
  labs(title = "Baseline: Poor vs Good Response",
       subtitle = "log2(OR) of positive cells. Red = Poor-enriched  Blue = Good-enriched.  *FDR<0.05  **FDR<0.01  ***FDR<0.001",
       x = NULL, y = NULL) +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7.5),
        axis.text.y = element_text(size = 7.5))
ggsave(file.path(FIG, "Fig6c_baseline_Tier123.pdf"), g_6c,
       width = 240, height = 120, units = "mm", limitsize = FALSE)

# (elided: Other-tier top-30 heatmap [sc_29 FigB -> Supp Fig S1] and four-group
#  % positive-cell panel [sc_29 FigC]; identical to Fig 6c construction.)

# -- 10.4 Fig 6d: paired delta panel ------------------------------------------
# (elided: drawn in Section 5.1 and saved as Fig6d_paired_delta.pdf; the
#  published panel is the same plot with per-pair mean +/- SE.)

# -- 10.5 Fig 6e: unified evidence matrix, 13 Tier 1-3 genes x 6 criteria -----
# Criteria (paper Results 2.5): E1 MR IVW FDR < 0.05; E2 coloc
# PP.H4/(PP.H3+PP.H4) > 0.7; E3 SMR (p < 0.05) + HEIDI (p > 0.05);
# E4 cross-database overlap (Soskic + DICE + eQTLGen); T1 tractability;
# T2 AI-screening feasibility (curated per the paper's evidence tables).
# TLR1 satisfies all 6 criteria (first); IL1R2 only E1 (last).
coloc_tab <- fread(file.path(RES, "06_coloc.tsv"))
coloc_sig_genes <- coloc_tab[coloc_sig == TRUE, unique(gene)]
smr_tab <- fread(file.path(RES, "A08_SMR_HEIDI.tsv"))
smr_pass_genes <- smr_tab[pass_both == TRUE, unique(gene)]
cross_tab <- fread(file.path(RES, "B03_crossDB_summary_AR.tsv"))
triple_genes <- cross_tab[n_db == 3, unique(gene)]
gm <- fread(file.path(RES, "gene_symbol_to_ensembl.tsv"))
setnames(gm, c("hgnc_symbol", "ensembl_gene_id"), c("gene", "ensg"))

evidence <- valid_genes[tier %in% c("Tier1", "Tier2", "Tier3"),
                        .(gene = gene_symbol, tier)]
evidence <- merge(evidence, gm, by = "gene", all.x = TRUE)
evidence[, E1 := 1L]                          # all Tier 1-3 genes: MR FDR < 0.05
evidence[, E2 := as.integer(ensg %in% coloc_sig_genes)]
evidence[, E3 := as.integer(ensg %in% smr_pass_genes)]
evidence[, E4 := as.integer(ensg %in% triple_genes)]
# T1/T2 curation follows the paper's evidence scoring tables
T1 <- c(TSLP = 0L, ORMDL3 = 0L, SLC25A46 = 1L, TLR1 = 1L, IL18R1 = 0L, IL18RAP = 0L,
        IKZF3 = 1L, CD247 = 1L, PDCD1 = 0L, D2HGDH = 0L, STAT6 = 1L, LMBRD2 = 0L, IL1R2 = 0L)
T2 <- c(TSLP = 1L, ORMDL3 = 0L, SLC25A46 = 0L, TLR1 = 1L, IL18R1 = 0L, IL18RAP = 0L,
        IKZF3 = 1L, CD247 = 1L, PDCD1 = 0L, D2HGDH = 0L, STAT6 = 0L, LMBRD2 = 0L, IL1R2 = 0L)
evidence[, `:=`(T1 = T1[gene], T2 = T2[gene])]
evidence[, score := E1 + E2 + E3 + E4 + T1 + T2]
setorder(evidence, -score, tier, gene)
evidence[, gene_ord := factor(gene, levels = rev(gene))]

dim_labs <- c("E1 MR", "E2 Coloc", "E3 SMR+HEIDI", "E4 Cross-DB", "T1", "T2")
crit_cols <- c("E1", "E2", "E3", "E4", "T1", "T2")
plot_df <- evidence[, .(gene_ord = rep(gene_ord, each = 6),
                        tier = rep(tier, each = 6),
                        score = rep(score, each = 6))]
plot_df[, dim := rep(dim_labs, times = nrow(evidence))]
plot_df[, value := unlist(evidence[, .SD, .SDcols = crit_cols])]
plot_df[, dim := factor(dim, levels = rev(dim_labs))]

g_6e <- ggplot(plot_df, aes(x = dim, y = gene_ord)) +
  geom_tile(aes(fill = factor(value)), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(value == 1, "+", "")), size = 3.5,
            colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("0" = "#E8E8E8", "1" = NPG_BLUE), guide = "none") +
  geom_text(data = unique(plot_df[, .(gene_ord, score)]),
            aes(x = 6.8, y = gene_ord, label = sprintf("%d/6", score)),
            size = 3.2, fontface = "bold", hjust = 0) +
  annotate("text", x = 6.8, y = nrow(evidence) + 0.5, label = "Score",
           size = 3, fontface = "bold", hjust = 0) +
  scale_x_discrete(position = "top", expand = expansion(add = c(0.5, 1.5))) +
  labs(title = "MR-scRNA unified evidence (13 genes x 6 criteria)",
       x = NULL, y = NULL) +
  theme_pub(9) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 9, face = "italic"),
        axis.ticks = element_blank(), plot.margin = margin(5, 40, 5, 30))
ggsave(file.path(FIG, "Fig6e_evidence_matrix.pdf"), g_6e,
       width = 200, height = 120, units = "mm", limitsize = FALSE)

# -- 10.6 Fig 6f: miloR beeswarm ----------------------------------------------
# (elided: drawn in Section 7 and saved as Fig6f_miloR_beeswarm.pdf.)

# -- 10.7 Supplementary figure panel map ---------------------------------------
# Supp S1: Other-tier top-30 baseline heatmap (10.3 pattern, tier == "Other")
# Supp S2: GSE273975 Tier 1-3 Pre vs Post Ro/e (sc_31 pattern)
# Supp S3: Cross-dataset direction scatter (sc_33 pattern)
# Supp S4: miloR neighborhood graph (sc_30 pattern)
# Supp S5: MR gene dotplot (Section 3) | Supp S6: pseudobulk volcano (Section 6)
# Supp S7: TCR figures (Section 8) | Supp S8: CellChat MR-LR bubble (Section 9.1)
# Supp S9: TF regulon dotplots (Section 9.2) | Supp S10: Mito-MR + PDCD1 (Section 9.3)

# =============================================================================
# End of Results 2.5 script
# =============================================================================
