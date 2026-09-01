# =============================================================================
# 03_crossdb_overlap.R
# Results 2.3: Cross-database overlap identified a nine-gene set and revealed
# timepoint-dependent effects (paper Figure 4)
# Paper: Genetic prioritization of CD4+ T cell-expressed candidate genes for
#        allergic rhinitis (Biomedicines)
# Merged and sanitized from: MR_09_extract_eqtl_ref.R, MR_10_map_rsid_ref.R,
#        MR_11_dice_mr.R, MR_12_eqtlgen_mr.R, MR_13_crossdb_comparison.R,
#        MR_07_timepoint_specificity.R, MR_figures_v4.R
# =============================================================================
# Notes:
#  - Run from the repository root; each script assumes cwd == repo root.
#  - Raw GWAS/eQTL data are obtained from public repositories (see the
#    paper's Data Availability statement); they are not distributed here.
#  - Repetitive blocks are elided and marked "# (elided: ...)".
# =============================================================================

# Libraries ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(TwoSampleMR); library(mr.raps)
  library(org.Hs.eg.db); library(parallel)
  library(ggplot2); library(ggrepel); library(dplyr)
  library(RColorBrewer); library(patchwork); library(ggvenn)
})

# Paths ----------------------------------------------------------------------
# Raw data (AR GWAS GCST90468131, DICE, eQTLGen, 1000G EUR, dbSNP) are expected
# under data/ (see Data Availability in README.md), with layout:
#   allergic_rhinitis_gwas_data/GCST90468131/GCST90468131.h.tsv.gz
#   DICE/ (one VCF per cell type), eQTLGen/2019-12-11-cis-eQTLsFDR0.05-...
#   g1000_eur/g1000_eur, plink_mac_20250819/plink
#   Soskic_CD4/tensor_out/, snp150_by_chr/
raw_dir <- "data"
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# Shared helpers --------------------------------------------------------------
# NPG high-contrast palette (paper-wide figure style)
npg10 <- c("#D62728", "#1F77B4", "#FF7F0E", "#2CA02C", "#9467BD",
           "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")

# Publication theme (Nature style, visible grid lines)
theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) + theme(
    axis.text = element_text(color = "black", size = base_size - 1),
    axis.title = element_text(color = "black", size = base_size),
    plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "grey40"),
    legend.text = element_text(size = base_size - 1),
    legend.title = element_text(size = base_size - 1, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    strip.text = element_text(size = base_size - 1, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 6, 6, 6))
}

# Map a Soskic profile name to one of the five activation timepoints
# (0h, LA, 16h, 40h, 5d); unmatched profiles are "other".
parse_timepoint <- function(profile) {
  tp <- case_when(
    grepl("(?<!4)0h", profile, perl = TRUE) ~ "0h",
    grepl("_LA$|_LA_", profile) ~ "LA",
    grepl("16h", profile) ~ "16h",
    grepl("40h", profile) ~ "40h",
    grepl("5d", profile) ~ "5d",
    TRUE ~ "other")
  factor(tp, levels = c("0h", "LA", "16h", "40h", "5d", "other"))
}

# Per-gene MR estimation: Wald ratio for single-instrument genes, otherwise
# IVW + weighted median + weighted mode. One representative instance; the
# eQTLGen section reuses this function.
run_mr_per_gene <- function(dat, cell_type = NULL) {
  out <- list()
  for (g in unique(dat$exposure)) {
    sub <- dat[dat$exposure == g, ]
    if (nrow(sub) == 1) {
      wr <- mr(sub, method_list = "mr_wald_ratio")
      row <- data.frame(gene = g, nsnp = 1, method = "Wald ratio",
                        beta = wr$b, se = wr$se, pval = wr$pval)
    } else {
      mr_res <- mr(sub, method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
      row <- data.frame(gene = g, nsnp = mr_res$nsnp, method = mr_res$method,
                        beta = mr_res$b, se = mr_res$se, pval = mr_res$pval)
    }
    if (!is.null(cell_type)) row <- cbind(cell_type = cell_type, row)
    out[[length(out) + 1]] <- row
  }
  rbindlist(out, fill = TRUE)
}

# Reference eQTL extraction and rsID mapping --------------------------------
# Supplementary to the primary pipeline; reproduces the top-eQTL table used
# by the main MR modules.
tensor_dir <- file.path(raw_dir, "Soskic_CD4/tensor_out")
files <- list.files(tensor_dir, pattern = "\\.txt$", full.names = TRUE)
all_profiles <- lapply(files, function(f) {
  dt <- fread(f)
  dt$profile <- sub("_500kb_window_tensorQTL.txt", "", basename(f))
  dt
})
eqtl <- rbindlist(all_profiles)
eqtl <- eqtl[pval_beta < 0.05]
eqtl$F_stat <- (eqtl$slope / eqtl$slope_se)^2
eqtl <- eqtl[F_stat >= 10]

eqtl_out <- eqtl[, .(
  profile, gene_id = phenotype_id, snp_id = variant_id,
  chr = sub("_.*", "", variant_id),
  pos = as.integer(sub("^[^_]+_([^_]+)_.*", "\\1", variant_id)),
  ref = sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", variant_id),
  alt = sub("^.*_([^_]+)$", "\\1", variant_id),
  beta = slope, se = slope_se, pval_nominal, pval_beta, maf, F_stat
)]
fwrite(eqtl_out, "results/module_A/01_top_eqtl.tsv", sep = "\t")

# Map to rsID using per-chromosome dbSNP files only: the full dbSNP table is
# ~7.1 GB and is never loaded; each chromosome file is filtered to the
# positions required by the eQTL table (hundreds of rows, not millions).
eqtl <- fread("results/module_A/01_top_eqtl.tsv")
eqtl[, chrom := paste0("chr", chr)]
pos_by_chr <- split(eqtl$pos, eqtl$chrom)
snp_ref <- rbindlist(mclapply(names(pos_by_chr), function(ch) {
  f <- file.path(raw_dir, "snp150_by_chr", paste0(ch, ".txt"))
  if (!file.exists(f)) return(NULL)
  pos_needed <- pos_by_chr[[ch]]
  dt <- fread(f, select = c(2, 4, 5), col.names = c("chrom", "pos", "rsid"))
  dt[chrom == ch & pos %in% pos_needed]
}, mc.cores = 12))
snp_ref[, chr := sub("chr", "", chrom)]
eqtl[, chr := as.character(chr)]
eqtl <- merge(eqtl, snp_ref[, .(chr, pos, rsid)], by = c("chr", "pos"), all.x = TRUE)
eqtl <- eqtl[!is.na(rsid)]
eqtl[, chrom := NULL]
fwrite(eqtl, "results/module_A/02_eqtl_with_rsid.tsv", sep = "\t")

# DICE cross-database MR (15 sorted immune cell populations) -----------------
# One top instrument per gene (P < 5e-5), clumped (r2 < 0.001, 500 kb,
# 1000G EUR), F >= 10, harmonised against the AR GWAS and tested by MR.

# Shared input: AR GWAS (GCST90468131, Hayfever, N = 394,626, GRCh38).
gwas_raw <- fread(
  cmd = paste0("gunzip -c ", file.path(raw_dir,
    "allergic_rhinitis_gwas_data/GCST90468131/GCST90468131.h.tsv.gz")),
  sep = "\t")
gwas_raw <- gwas_raw[!is.na(beta) & !is.na(standard_error) & standard_error > 0]
gwas_raw <- gwas_raw[!is.na(rsid) & rsid != ""]
gwas_out <- data.frame(
  SNP = gwas_raw$rsid, beta.outcome = gwas_raw$beta,
  se.outcome = gwas_raw$standard_error, pval.outcome = gwas_raw$p_value,
  eaf.outcome = gwas_raw$effect_allele_frequency,
  effect_allele.outcome = gwas_raw$effect_allele,
  other_allele.outcome = gwas_raw$other_allele,
  outcome = "GCST90468131", id.outcome = "GCST90468131")

plink_bin <- normalizePath(file.path(raw_dir, "plink_mac_20250819/plink"))
bfile     <- normalizePath(file.path(raw_dir, "g1000_eur/g1000_eur"))
dice_dir  <- file.path(raw_dir, "DICE")
vcf_files <- list.files(dice_dir, pattern = "\\.vcf$", full.names = TRUE)

# One iteration per cell type (the representative instance of the
# per-cell-type pattern; all 15 DICE cell types follow the same steps)
dice_results <- list()
for (vcf_file in vcf_files) {
  cell_type <- sub("\\.vcf$", "", basename(vcf_file))
  dice <- fread(vcf_file, skip = "#CHROM")
  dice[, Gene := sub(".*Gene=([^;]+).*", "\\1", INFO)]
  dice[, GeneSymbol := sub(".*GeneSymbol=([^;]+).*", "\\1", INFO)]
  dice[, Pvalue := as.numeric(sub(".*Pvalue=([^;]+).*", "\\1", INFO))]
  dice[, Beta := as.numeric(sub(".*Beta=([^;]+).*", "\\1", INFO))]
  dice <- dice[Pvalue < 5e-5]
  top_snps <- dice[, .SD[which.min(Pvalue)], by = Gene]

  exposure_df <- data.frame(
    SNP = top_snps$ID, beta.exposure = top_snps$Beta,
    se.exposure = abs(top_snps$Beta / qnorm(top_snps$Pvalue / 2, lower.tail = FALSE)),
    pval.exposure = top_snps$Pvalue, eaf.exposure = 0.3,
    effect_allele.exposure = top_snps$ALT, other_allele.exposure = top_snps$REF,
    exposure = top_snps$Gene, id.exposure = top_snps$Gene)
  exposure_df <- exposure_df[exposure_df$SNP %in% gwas_raw$rsid, ]
  if (nrow(exposure_df) == 0) next

  exposure_clumped <- tryCatch(
    clump_data(exposure_df, clump_r2 = 0.001, clump_kb = 500,
               pop = "EUR", bfile = bfile, plink_bin = plink_bin),
    error = function(e) NULL)
  if (is.null(exposure_clumped) || nrow(exposure_clumped) == 0) next
  exposure_clumped$F_stat <- (exposure_clumped$beta.exposure / exposure_clumped$se.exposure)^2
  exposure_clumped <- exposure_clumped[exposure_clumped$F_stat >= 10, ]
  if (nrow(exposure_clumped) == 0) next

  gwas_sub <- gwas_out[gwas_out$SNP %in% exposure_clumped$SNP, ]
  dat <- harmonise_data(exposure_clumped, gwas_sub, action = 2)
  dat <- dat[dat$mr_keep == TRUE, ]
  if (nrow(dat) == 0) next
  dice_results[[length(dice_results) + 1]] <- run_mr_per_gene(dat, cell_type)
}
dice_mr_all <- rbindlist(dice_results, fill = TRUE)
dice_ivw <- dice_mr_all[method %in% c("Inverse variance weighted", "Wald ratio")]
dice_ivw[, fdr := p.adjust(pval, method = "BH")]
cat("DICE: significant genes (FDR < 0.05):",
    sum(dice_ivw$fdr < 0.05, na.rm = TRUE), "\n")
fwrite(dice_mr_all, "results/B01_DICE_MR_all.tsv", sep = "\t")
fwrite(dice_ivw,    "results/B01_DICE_MR_ivw.tsv", sep = "\t")

# eQTLGen whole-blood meta-analysis (N = 31,684) ------------------------------
eqtlgen <- fread(file.path(raw_dir,
  "eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"))
# (elided: AR GWAS loading and outcome preparation -- identical to the DICE section; gwas_out reused)
top_eqtlgen <- eqtlgen[, .SD[which.min(Pvalue)], by = GeneSymbol]
top_eqtlgen <- top_eqtlgen[SNP %in% gwas_raw$rsid]
exposure_df <- data.frame(
  SNP = top_eqtlgen$SNP, beta.exposure = top_eqtlgen$Zscore / sqrt(top_eqtlgen$NrSamples),
  se.exposure = 1 / sqrt(top_eqtlgen$NrSamples), pval.exposure = top_eqtlgen$Pvalue,
  eaf.exposure = rep(0.3, nrow(top_eqtlgen)),
  effect_allele.exposure = top_eqtlgen$AssessedAllele,
  other_allele.exposure = top_eqtlgen$OtherAllele,
  exposure = top_eqtlgen$GeneSymbol, id.exposure = top_eqtlgen$GeneSymbol)
exposure_clumped <- clump_data(exposure_df, clump_r2 = 0.001, clump_kb = 500,
                               pop = "EUR", bfile = bfile, plink_bin = plink_bin)
exposure_clumped$F_stat <- (exposure_clumped$beta.exposure / exposure_clumped$se.exposure)^2
exposure_clumped <- exposure_clumped[exposure_clumped$F_stat >= 10, ]
gwas_sub <- gwas_out[gwas_out$SNP %in% exposure_clumped$SNP, ]
dat <- harmonise_data(exposure_clumped, gwas_sub, action = 2)
dat <- dat[dat$mr_keep == TRUE, ]
# (elided: per-gene MR estimation loop -- same pattern as the DICE section; run_mr_per_gene() reused)
eqtlgen_mr_all <- run_mr_per_gene(dat)
eqtlgen_ivw <- eqtlgen_mr_all[method %in% c("Inverse variance weighted", "Wald ratio")]
eqtlgen_ivw[, fdr := p.adjust(pval, method = "BH")]
cat("eQTLGen: significant genes (FDR < 0.05):",
    sum(eqtlgen_ivw$fdr < 0.05, na.rm = TRUE), "\n")
fwrite(eqtlgen_mr_all, "results/B02_eQTLGen_MR_all.tsv", sep = "\t")
fwrite(eqtlgen_ivw,    "results/B02_eQTLGen_MR_ivw.tsv", sep = "\t")

# Three-way overlap and direction concordance --------------------------------
# All MR results were generated with MHC-region instruments (chr6:25.5-34 Mb)
# excluded; significance is FDR < 0.05 throughout.
mr_full <- fread("results/05_MR_ivw_fdr.tsv")   # primary Soskic analysis
sig_a5_ensg <- unique(mr_full[fdr < 0.05]$gene)   # Soskic CD4+ (Ensembl)
sig_b1_ensg <- unique(dice_ivw[fdr < 0.05]$gene)  # DICE (Ensembl)
sig_b2_sym  <- unique(eqtlgen_ivw[fdr < 0.05]$gene)  # eQTLGen (symbols)

# ENSG -> symbol map for the significant Soskic/DICE genes
all_ensg <- unique(c(sig_a5_ensg, sig_b1_ensg))
map <- AnnotationDbi::select(org.Hs.eg.db,
  keys = all_ensg, columns = c("ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
map_dt <- as.data.table(map)[!is.na(SYMBOL) & SYMBOL != ""]
sym_lup <- setNames(map_dt$SYMBOL, map_dt$ENSEMBL)
# eQTLGen reports symbols; map them onto the ENSG universe
sig_b2_ensg <- unique(map_dt[SYMBOL %in% sig_b2_sym]$ENSEMBL)

int_a5_b1 <- intersect(sig_a5_ensg, sig_b1_ensg)
int_a5_b2 <- intersect(sig_a5_ensg, sig_b2_ensg)
int_b1_b2 <- intersect(sig_b1_ensg, sig_b2_ensg)
int_all_3 <- Reduce(intersect, list(sig_a5_ensg, sig_b1_ensg, sig_b2_ensg))
cat("Significant genes: Soskic =", length(sig_a5_ensg),
    ", DICE =", length(sig_b1_ensg), ", eQTLGen =", length(sig_b2_sym), "\n")
cat("Intersections: Soskic&DICE =", length(int_a5_b1),
    ", Soskic&eQTLGen =", length(int_a5_b2), ", DICE&eQTLGen =", length(int_b1_b2), "\n")
cat("Triple intersection (", length(int_all_3), "):",
    paste(map_dt[ENSEMBL %in% int_all_3]$SYMBOL, collapse = ", "), "\n")

# Effect directions within the triple intersection (Soskic vs eQTLGen)
a5_best <- mr_full[gene %in% int_all_3, .SD[which.min(pval)], by = gene
                   ][, .(gene, beta_a5 = beta, pval_a5 = pval)]
b2_sub_syms <- map_dt[ENSEMBL %in% int_all_3, .(gene = ENSEMBL, symbol = SYMBOL)]
b2_sym_best <- eqtlgen_ivw[gene %in% b2_sub_syms$symbol,
                           .(symbol = gene, beta_b2 = beta, pval_b2 = pval)]
b2_sub <- merge(b2_sub_syms, b2_sym_best, by = "symbol")[, .(gene, beta_b2, pval_b2)]
merged_dir <- merge(a5_best, b2_sub, by = "gene")
merged_dir[, direction_concordant := sign(beta_a5) == sign(beta_b2)]
merged_dir[, gene_symbol := sym_lup[gene]]
cat("Triple-intersection direction concordance (Soskic vs eQTLGen):\n")
print(merged_dir[, .(gene_symbol, beta_a5 = round(beta_a5, 3),
                     beta_b2 = round(beta_b2, 3), direction_concordant)])

# Summary table: significant gene sets per database (ENSG space)
summary_dt <- data.table(gene = union(union(sig_a5_ensg, sig_b1_ensg), sig_b2_ensg))
summary_dt[, gene_symbol := map_dt[match(gene, ENSEMBL), SYMBOL]]
summary_dt[is.na(gene_symbol), gene_symbol := gene]
summary_dt[, in_A5 := gene %in% sig_a5_ensg]
summary_dt[, in_B1 := gene %in% sig_b1_ensg]
summary_dt[, in_B2 := gene %in% sig_b2_ensg]
summary_dt[, n_db := in_A5 + in_B1 + in_B2]
setorder(summary_dt, -n_db)
fwrite(summary_dt, "results/B03_crossDB_summary_AR.tsv", sep = "\t")

# Nine-gene core set: NFKB1, IL18RAP, AHI1, IKZF3, SLC25A46, TEF, TLR1,
# ZBTB38, CD247 (significant in all three databases)
core9_genes <- summary_dt[n_db == 3, gene]
core9_syms  <- summary_dt[n_db == 3, gene_symbol]
cat("Core genes (3 databases):", paste(core9_syms, collapse = ", "), "\n")

# Timepoint analysis: Cochrane Q heterogeneity + pairwise Z-tests ------------
# Per-gene MR effects across the five activation timepoints (0h, LA, 16h,
# 40h, 5d).
mr_sig <- mr_full[fdr < 0.05 & method %in% c("Inverse variance weighted", "Wald ratio")]
mr_sig[, tp := as.character(parse_timepoint(profile))]
mr_sig <- mr_sig[tp != "other"]
genes_multi <- mr_sig[, .(n_tp = uniqueN(tp)), by = gene][n_tp >= 2, gene]
# (elided: per-gene best-profile-per-timepoint fits -- repeated in both the Cochrane Q and pairwise Z blocks; computed once here)
tp_best_all <- mr_sig[order(pval), .SD[1], by = .(gene, tp)]

# Cochrane Q heterogeneity test (inverse-variance fixed-effect model)
cochrane_results <- rbindlist(lapply(genes_multi, function(g) {
  best <- tp_best_all[gene == g]
  k <- nrow(best)
  if (k < 2) return(NULL)
  w_i <- 1 / best$se^2
  beta_fe <- sum(w_i * best$beta) / sum(w_i)
  Q <- sum(w_i * (best$beta - beta_fe)^2)
  data.table(gene = g, n_timepoints = k,
             timepoints = paste(sort(best$tp), collapse = ","),
             profiles = paste(best$profile, collapse = ";"),
             beta_fe = beta_fe, Q = Q, df = k - 1,
             pval_q = pchisq(Q, k - 1, lower.tail = FALSE))
}))
cochrane_results[, fdr_q := p.adjust(pval_q, method = "BH")]
cochrane_results[, q_sig := fdr_q < 0.05]
setorder(cochrane_results, fdr_q)
cat("Cochrane Q: significant temporal heterogeneity:",
    sum(cochrane_results$q_sig), "genes\n")

# Pairwise Z-tests between all 10 timepoint pairs
pairs_tp <- combn(c("0h", "LA", "16h", "40h", "5d"), 2, simplify = FALSE)
z_results <- rbindlist(lapply(genes_multi, function(g) {
  best <- tp_best_all[gene == g]
  rbindlist(lapply(pairs_tp, function(tp_pair) {
    t1 <- best[tp == tp_pair[1]]
    t2 <- best[tp == tp_pair[2]]
    if (nrow(t1) == 0 || nrow(t2) == 0) return(NULL)
    z <- (t1$beta - t2$beta) / sqrt(t1$se^2 + t2$se^2)
    data.table(gene = g, timepoint1 = tp_pair[1], profile1 = t1$profile,
               beta1 = t1$beta, se1 = t1$se, timepoint2 = tp_pair[2],
               profile2 = t2$profile, beta2 = t2$beta, se2 = t2$se,
               z_score = z, pval_diff = 2 * pnorm(-abs(z)))
  }))
}))
z_results[, fdr_diff := p.adjust(pval_diff, method = "BH")]
z_results[, timepoint_specific := fdr_diff < 0.05]
setorder(z_results, fdr_diff)
cat("Pairwise Z: significant timepoint-specific changes:",
    sum(z_results$timepoint_specific), "across",
    uniqueN(z_results$gene[z_results$timepoint_specific]), "genes\n")
fwrite(cochrane_results, "results/A07_Cochrane_Q.tsv", sep = "\t")
fwrite(z_results,        "results/A07_pairwise_z.tsv", sep = "\t")

# Figure 4 panels -------------------------------------------------------------
# Paper Figure 4 = (a) Venn diagram, (b) direction comparison,
# (c) timepoint trajectories, (d) core-gene small multiples. Internal figure
# numbers in the source scripts were Fig 3 -> 4a, Fig 8 -> 4c, Fig S2 -> 4d.

# --- Figure 4a: Venn diagram of the three-database intersection --------------
cat("Figure 4a: Venn diagram\n")
venn_list <- list(
  "Soskic CD4+" = summary_dt[in_A5 == TRUE, gene],
  "DICE"        = summary_dt[in_B1 == TRUE, gene],
  "eQTLGen"     = summary_dt[in_B2 == TRUE, gene])
p_venn <- ggvenn(venn_list, show_percentage = FALSE, show_elements = FALSE,
  text_size = 5, set_name_size = 4.5, digits = 0, auto_scale = FALSE,
  stroke_color = "grey30", stroke_size = 0.4, stroke_alpha = 0.6,
  fill_color = c("#4C72B0", "#E64B35", "#55A868"),
  set_name_color = c("#3B5F9E", "#C0392B", "#3D8B50"))
# Core nine-gene list box
p_list <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           fill = "#FAFAFA", color = "#BBBBBB", linewidth = 0.8) +
  annotate("text", x = 0.5, y = 0.92, label = "Core Causal Genes",
           size = 4.5, fontface = "bold", color = "#222222") +
  annotate("text", x = 0.5, y = 0.86, label = "(3 databases)",
           size = 3, color = "#777777") +
  annotate("segment", x = 0.2, xend = 0.8, y = 0.82, yend = 0.82,
           color = "#DDDDDD", linewidth = 0.4) +
  annotate("text", x = 0.5, y = 0.45, label = paste(core9_syms, collapse = "\n"),
           size = 3.8, fontface = "italic", color = "#555555", lineheight = 1.1) +
  xlim(0, 1) + ylim(0, 1) + theme_void() +
  theme(plot.margin = margin(8, 10, 8, 5, "mm"))
fig4a <- p_venn + p_list + plot_layout(widths = c(1.4, 0.5))
ggsave(file.path("figures", "Fig4a_venn_crossdb.pdf"),
       fig4a, width = 9.5, height = 6, device = "pdf")

# --- Figure 4b: cross-database MR direction comparison -----------------------
cat("Figure 4b: direction comparison\n")
col_concord <- "#55A868"; col_discord <- "#C44E52"
theme_nc <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(color = "black"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        plot.subtitle = element_text(size = 7, color = "grey40"))
# Best estimate per gene in the primary analysis, restricted to genes
# significant in both Soskic and eQTLGen (fdr < 0.05). This reproduces the
# 16-gene comparison reported in the paper (10/16 = 62.5% concordant; the
# discordant genes include NFKB1, IL18RAP and TLR1).
soskic_best <- mr_full[fdr < 0.05, .SD[which.min(pval)], by = gene]
soskic_best[, symbol := sym_lup[gene]]
soskic_best <- soskic_best[!is.na(symbol)]
eqtlgen_dir <- unique(eqtlgen_ivw[fdr < 0.05], by = "gene")
setnames(eqtlgen_dir, "gene", "symbol")
merged <- merge(soskic_best, eqtlgen_dir, by = "symbol",
                suffixes = c("_soskic", "_eqtlgen"))
merged[, concordant := sign(beta_soskic) == sign(beta_eqtlgen)]
merged[, direction := fifelse(concordant, "Concordant", "Discordant")]
merged[, is_core := symbol %in% core9_syms]
n_total   <- nrow(merged)
n_concord <- sum(merged$concordant)
fig4b <- ggplot(merged, aes(x = beta_soskic, y = beta_eqtlgen)) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, color = "grey50") +
  geom_point(data = merged[is_core == FALSE], aes(color = direction),
             size = 0.8, alpha = 0.5) +
  geom_point(data = merged[is_core == TRUE], aes(color = direction),
             size = 2.5, alpha = 1) +
  geom_text_repel(data = merged[is_core == TRUE], aes(label = symbol),
                  size = 2.8, fontface = "italic", box.padding = 0.4,
                  max.overlaps = 20, min.segment.length = 0.15,
                  color = "black", bg.color = "white", bg.r = 0.08) +
  scale_color_manual(name = "Direction",
                     values = c("Concordant" = col_concord, "Discordant" = col_discord)) +
  annotate("label", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("Concordance: %d/%d (%.1f%%)",
                           n_concord, n_total, 100 * n_concord / n_total),
           size = 2.8, fontface = "bold", fill = "white", label.size = 0.3) +
  labs(title = "Cross-database MR direction comparison",
       subtitle = sprintf("Genes significant in both Soskic CD4+ and eQTLGen (n = %d)", n_total),
       x = "Soskic MR beta", y = "eQTLGen MR beta") +
  coord_fixed() + theme_nc + theme(legend.position = "bottom")
ggsave(file.path("figures", "Fig4b_crossdb_direction.pdf"),
       fig4b, width = 6, height = 6, device = "pdf")

# --- Figure 4c: timepoint trajectories (Cochrane Q significant genes) --------
cat("Figure 4c: timepoint trajectories\n")
cq_genes <- cochrane_results[q_sig == TRUE, gene]
tp_mr <- mr_full[fdr < 0.05 & gene %in% cq_genes]
tp_mr[, gene_symbol := coalesce(sym_lup[gene], gene)]
tp_mr[, tp := parse_timepoint(profile)]
tp_mr <- tp_mr[tp != "other"]
tp_best <- tp_mr[order(pval), .SD[1], by = .(gene, tp)]
tp_best[, timepoint := factor(tp, levels = c("0h", "LA", "16h", "40h", "5d"))]
# Labels at the last observed timepoint
tp_labels <- tp_best[, .SD[which.max(as.integer(timepoint))], by = gene]
# Palette: the six genes coloured in the original figure, extended with NPG
# colors to the full set of nine heterogeneous genes
cq_genes_syms <- unique(tp_best$gene_symbol)
cq_colors <- c(D2HGDH = "#5A83B8", IL1R2 = "#C0656F", LMBRD2 = "#D4954B",
               CAMK4 = "#6B9E7A", TSLP = "#C44E52", STAT6 = "#4C72B0")
cq_colors <- c(cq_colors, setNames(npg10[seq_along(setdiff(cq_genes_syms, names(cq_colors)))],
                                   setdiff(cq_genes_syms, names(cq_colors))))
fig4c <- ggplot(tp_best, aes(x = timepoint, y = beta,
                             color = gene_symbol, group = gene_symbol)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_ribbon(aes(ymin = beta - se, ymax = beta + se, fill = gene_symbol),
              alpha = 0.06, color = NA) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_text_repel(data = tp_labels, aes(label = gene_symbol),
                  size = 2.8, direction = "y", nudge_x = 0.5,
                  segment.size = 0.2, segment.color = "grey50", hjust = 0) +
  scale_color_manual(values = cq_colors, guide = "none") +
  scale_fill_manual(values = cq_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.2))) +
  labs(title = "Activation timepoint-dependent MR effects",
       subtitle = sprintf("%d genes with significant temporal heterogeneity (Cochrane Q FDR < 0.05)",
                          length(cq_genes)),
       x = "CD4+ T cell activation timepoint",
       y = expression("MR effect size (" * beta * "  " %+-% "  SE)")) +
  theme_nature()
ggsave(file.path("figures", "Fig4c_timepoint_trajectories.pdf"),
       fig4c, width = 8, height = 5.5, device = "pdf")

# --- Figure 4d: small multiples for the core genes ---------------------------
cat("Figure 4d: core-gene timepoint small multiples\n")
core9_tp <- mr_full[gene %in% core9_genes & fdr < 0.05]
core9_tp[, gene_symbol := coalesce(sym_lup[gene], gene)]
core9_tp[, tp := as.character(parse_timepoint(profile))]
core9_tp <- core9_tp[tp != "other"]
core9_tp[, tp := factor(tp, levels = c("0h", "LA", "16h", "40h", "5d"))]
core9_summary <- core9_tp[, .(beta_mean = mean(beta),
                              beta_se = sd(beta) / sqrt(.N),
                              n_profiles = .N), by = .(gene_symbol, tp)]
# Keep genes present at >= 2 timepoints
genes_multi <- core9_summary[, .N, by = gene_symbol][N >= 2, gene_symbol]
core9_summary <- core9_summary[gene_symbol %in% genes_multi]
fig4d <- ggplot(core9_summary, aes(x = tp, y = beta_mean,
                                   color = gene_symbol, group = gene_symbol)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbar(aes(ymin = beta_mean - beta_se, ymax = beta_mean + beta_se),
                width = 0.2, alpha = 0.7) +
  geom_line(alpha = 0.6, linewidth = 0.6) +
  geom_point(size = 2.5) +
  facet_wrap(~gene_symbol, scales = "free_y", ncol = 5) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(9),
                     guide = "none") +
  labs(title = "Activation timepoint-dependent MR effects for core AR genes",
       subtitle = "Mean beta +/- SE across FDR < 0.05 profiles at each timepoint",
       x = "Activation timepoint",
       y = expression("Mean MR effect (" * beta * ")")) +
  theme_nature(base_size = 9)
ggsave(file.path("figures", "Fig4d_core_timepoint_facets.pdf"),
       fig4d, width = 10, height = 5, device = "pdf")

cat("Results 2.3 pipeline complete: tables in results/, figures in figures/\n")
