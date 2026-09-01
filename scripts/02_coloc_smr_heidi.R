# =============================================================================
# 02_coloc_smr_heidi.R
# Results 2.2: Colocalization and SMR/HEIDI reduced the candidate set to a
# more credible subset (paper Figure 3)
# Paper: Genetic prioritization of CD4+ T cell-expressed candidate genes for
#        allergic rhinitis (Biomedicines)
# Merged and sanitized from: MR_06_coloc.R, MR_08_smr_heidi.R, MR_figures_v4.R
# =============================================================================
# Notes:
#  - Run from the repository root; each script assumes cwd == repo root.
#  - Raw GWAS/eQTL data are obtained from public repositories (see the
#    paper's Data Availability statement); they are not distributed here.
#  - Repetitive blocks are elided and marked "# (elided: ...)".
# =============================================================================

# ---- Libraries ---------------------------------------------------------------
library(data.table)
library(coloc)
library(arrow)
library(parallel)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(forcats)
})

# ---- Paths -------------------------------------------------------------------
# Raw data come from public repositories (GWAS Catalog GCST90468131, Soskic
# CD4+ T cell eQTL parquets, UCSC snp150 per-chromosome files) and are
# expected under data/ (see Data Availability in README.md).
RAW         <- "data"
gwas_path   <- file.path(RAW, "allergic_rhinitis_gwas_data/GCST90468131/GCST90468131.h.tsv.gz")
parquet_dir <- file.path(RAW, "Soskic_CD4")
split_dir   <- file.path(RAW, "snp150_by_chr")
RES  <- "results"
FIG  <- "figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
n_cores <- 4

# Disable data.table OpenMP threads before any fork to prevent OpenMP state
# corruption in forked child processes (mclapply uses fork on macOS/Linux).
setDTthreads(1)

# ---- Plot theme and NPG colors ------------------------------------------------
theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) + theme(
    axis.text        = element_text(color = "black", size = base_size - 1),
    axis.title       = element_text(color = "black", size = base_size),
    plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = base_size - 1, hjust = 0.5, color = "grey40"),
    legend.text      = element_text(size = base_size - 1),
    legend.title     = element_text(size = base_size - 1, face = "bold"),
    legend.key.size  = unit(0.4, "cm"),
    strip.text       = element_text(size = base_size - 1, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(6, 6, 6, 6))
}
col_risk    <- "#D62728"   # NPG red (risk-increasing)
col_protect <- "#1F77B4"   # NPG blue (protective)
col_core    <- "#FF7F0E"   # NPG orange (core gene highlight)
col_sig     <- "#2CA02C"   # NPG green (SMR + HEIDI pass)

# Map the 46 Soskic expression profiles to the five canonical activation time
# points: 0h (resting), LA (low activity), 16h, 40h, 5d.
parse_timepoint <- function(profile) {
  tp <- case_when(
    grepl("(?<!4)0h", profile, perl = TRUE) ~ "0h",
    grepl("_LA$|_LA_", profile)             ~ "LA",
    grepl("16h",  profile)                  ~ "16h",
    grepl("40h",  profile)                  ~ "40h",
    grepl("5d",   profile)                  ~ "5d",
    TRUE                                    ~ "other")
  factor(tp, levels = c("0h", "LA", "16h", "40h", "5d", "other"))
}

# ---- Shared input data --------------------------------------------------------
mr_full   <- fread(file.path(RES, "05_MR_ivw_fdr.tsv"))     # MR results
sig_genes <- unique(mr_full[fdr < 0.05, gene])              # 129 genes
eqtl_full <- fread(file.path(RES, "02_eqtl_with_rsid.tsv"),
                   select = c("gene_id", "chr", "pos"))     # eQTL gene positions
eqtl_full[, chr := as.character(chr)]
sym_dt  <- fread(file.path(RES, "gene_symbol_full.tsv"))    # Ensembl -> symbol
setnames(sym_dt, c("gene", "gene_symbol"))
sym_lup <- setNames(sym_dt$gene_symbol, sym_dt$gene)
# Nine core genes replicated in all three eQTL databases (Soskic, DICE,
# eQTLGen); taken from the cross-database replication results.
core9_genes <- fread(file.path(RES, "B03_crossDB_summary_AR.tsv"))[n_db == 3, gene]

# AR GWAS (GCST90468131, N = 394,626, ~20% cases): load the columns needed by
# colocalization (rsid, beta, se, eaf) and by SMR (chr, pos, alleles) in one
# parse-time column-selected read.
gwas_raw <- fread(
  cmd    = paste0("gunzip -c ", gwas_path),
  sep    = "\t",
  select = c("rsid", "chromosome", "base_pair_location", "effect_allele",
             "other_allele", "beta", "standard_error", "effect_allele_frequency"))
gwas_raw <- gwas_raw[!is.na(beta) & !is.na(standard_error) & standard_error > 0]
gwas_raw <- gwas_raw[!is.na(rsid) & rsid != ""]
gwas_raw[, varbeta := standard_error^2]
setkey(gwas_raw, rsid)

# =============================================================================
# 1. Bayesian colocalization of eQTL and AR GWAS signals
# =============================================================================
# For each MR-significant gene-profile pair, coloc.abf computes the posterior
# probability that both traits share a single causal variant (H4). We require
# the conditional ratio PP.H4/(PP.H3 + PP.H4) > 0.7 (reference-paper standard,
# NOT the raw PP.H4). In the current data 418 of 944 pairs (111 unique genes)
# meet this threshold.

# Gene-profile pairs to test, with the chromosome of each gene.
dat_all <- fread(file.path(RES, "04_clumped.tsv"),
                 select = c("exposure", "profile"))[exposure %in% sig_genes]
gene_profile_pairs <- unique(dat_all[, .(gene = exposure, profile)])
gene_chr <- unique(eqtl_full[gene_id %in% sig_genes, .(gene = gene_id, chr)])
gene_profile_pairs <- merge(gene_profile_pairs, gene_chr, by = "gene", all.x = TRUE)
gene_profile_pairs <- gene_profile_pairs[!is.na(chr)]

# Per-chromosome rsID lookup from snp150, pre-filtered with merged gene-region
# windows (+/- 500 kb) so the awk condition stays short in gene-dense regions.
chrs <- sort(unique(as.integer(gene_profile_pairs$chr)))
chrs <- chrs[!is.na(chrs)]
snp_by_chr   <- list()
region_rsids <- character(0)
for (chr_val in chrs) {
  snp_file <- file.path(split_dir, paste0("chr", chr_val, ".txt"))
  if (!file.exists(snp_file)) next
  gene_pos <- eqtl_full[chr == as.character(chr_val) & gene_id %in% sig_genes,
                        .(gene_id, pos)]
  if (nrow(gene_pos) == 0) next
  # Per-gene windows (+/- 500 kb), then merge overlapping intervals.
  wins <- gene_pos[, .(lo = max(1L, min(pos) - 500000L),
                        hi = max(pos) + 500000L), by = gene_id]
  wins <- wins[order(lo)]
  merged_wins <- wins[, {
    lo2 <- lo[1]; hi2 <- hi[1]; out <- list()
    for (k in seq_len(.N)) {
      if (lo[k] <= hi2) hi2 <- max(hi2, hi[k]) else {
        out[[length(out) + 1]] <- data.table(lo = lo2, hi = hi2)
        lo2 <- lo[k]; hi2 <- hi[k]
      }
    }
    rbindlist(c(out, list(data.table(lo = lo2, hi = hi2))))
  }]
  awk_cond <- paste(sprintf("($4 >= %d && $4 <= %d)", merged_wins$lo, merged_wins$hi),
                    collapse = " || ")
  snp_sub <- fread(
    cmd          = sprintf("awk -F'\t' '{if (%s) print}' '%s'", awk_cond, snp_file),
    select       = c(4, 5),
    col.names    = c("pos", "rsid"),
    showProgress = FALSE)
  setkey(snp_sub, pos)
  snp_by_chr[[as.character(chr_val)]] <- snp_sub
  region_rsids <- c(region_rsids, snp_sub$rsid)
}

# Subset the GWAS to region SNPs only (major memory reduction).
gwas_region <- gwas_raw[J(unique(region_rsids)), nomatch = 0]
setkey(gwas_region, rsid)   # J() subset drops the key; re-assert it
rm(region_rsids)
gc(verbose = FALSE)

# Parallelize by PROFILE: each forked worker loads its parquet exactly once;
# gwas_region and snp_by_chr are read-only in workers (copy-on-write sharing).
profiles <- unique(gene_profile_pairs$profile)
process_profile <- function(prof) {
  setDTthreads(1)   # re-assert in forked worker
  pf <- file.path(parquet_dir, paste0(prof, "_500kb_combined.parquet"))
  if (!file.exists(pf)) return(NULL)
  pq <- tryCatch({
    dt <- as.data.table(read_parquet(
      pf,
      col_select = c("phenotype_id", "variant_id", "slope", "slope_se", "maf")))
    dt[, chr_col := sub("_.*", "", variant_id)]
    dt[, pos     := as.integer(sub("^[^_]+_([^_]+)_.*", "\\1", variant_id))]
    setkey(dt, phenotype_id, chr_col)
    dt
  }, error = function(e) NULL)
  if (is.null(pq)) return(NULL)
  pairs_prof <- gene_profile_pairs[profile == prof]
  prof_results <- list()
  # Per-gene loop over the pairs of this profile (944 pairs overall in the
  # current data); each pair is processed identically.
  for (chr_str in unique(pairs_prof$chr)) {
    snp_ref <- snp_by_chr[[chr_str]]
    if (is.null(snp_ref) || nrow(snp_ref) == 0) next
    for (gene in pairs_prof[chr == chr_str, gene]) {
      eqtl_gene <- pq[.(gene, chr_str)]
      if (nrow(eqtl_gene) == 0) next
      region_start <- max(1L, min(eqtl_gene$pos) - 50000L)
      region_end   <- max(eqtl_gene$pos) + 50000L
      snp_region <- snp_ref[pos >= region_start & pos <= region_end]
      eqtl_gene  <- merge(eqtl_gene, snp_region[, .(pos, rsid)], by = "pos", all.x = TRUE)
      eqtl_gene  <- eqtl_gene[!is.na(rsid) & !is.na(slope) & !is.na(slope_se)]
      if (nrow(eqtl_gene) == 0) next
      gwas_sub    <- gwas_region[J(eqtl_gene$rsid), nomatch = 0]
      common_snps <- intersect(eqtl_gene$rsid, gwas_sub$rsid)
      if (length(common_snps) < 5) next
      eqtl_sub <- eqtl_gene[rsid %in% common_snps]
      gwas_sub <- gwas_sub[rsid %in% common_snps]
      coloc_res <- tryCatch(
        suppressMessages(coloc.abf(
          # eQTL exposure: quantitative trait, N = 119 (Soskic cohort).
          list(beta = eqtl_sub$slope, varbeta = eqtl_sub$slope_se^2,
               snp = eqtl_sub$rsid, type = "quant", N = 119L, MAF = eqtl_sub$maf),
          # AR GWAS outcome: case-control, N = 394,626, case fraction 0.20.
          list(beta = gwas_sub$beta, varbeta = gwas_sub$varbeta,
               snp = gwas_sub$rsid, type = "cc", N = 394626L, s = 0.20)))$summary,
        error = function(e) NULL)
      if (is.null(coloc_res)) next
      # unname() prevents the named vector from overriding data.table() names.
      prof_results[[length(prof_results) + 1]] <- data.table(
        gene = gene, profile = prof,
        PP.H0 = unname(coloc_res["PP.H0.abf"]), PP.H1 = unname(coloc_res["PP.H1.abf"]),
        PP.H2 = unname(coloc_res["PP.H2.abf"]), PP.H3 = unname(coloc_res["PP.H3.abf"]),
        PP.H4 = unname(coloc_res["PP.H4.abf"]))
    }
  }
  rm(pq)
  gc(verbose = FALSE)
  if (length(prof_results) == 0) return(NULL)
  rbindlist(prof_results)
}
results_list <- mclapply(profiles, process_profile, mc.cores = n_cores)

coloc_out <- rbindlist(Filter(function(x) is.data.table(x) || is.data.frame(x),
                              results_list))
if (nrow(coloc_out) == 0) stop("No colocalization results produced.")

# Conditional-ratio significance threshold (NOT raw PP.H4).
coloc_out[, coloc_sig   := PP.H4 / (PP.H3 + PP.H4) > 0.7]
coloc_out[, coloc_ratio := PP.H4 / (PP.H3 + PP.H4)]
n_coloc <- sum(coloc_out$coloc_sig)   # 418 of 944 pairs in the current data
cat("Coloc-significant pairs:", n_coloc, "of", nrow(coloc_out), "\n")
fwrite(coloc_out, file.path(RES, "06_coloc.tsv"), sep = "\t")

# =============================================================================
# 2. Colocalization support across the five activation timepoints
# =============================================================================
# For each of the 9 core genes, take the maximum PP.H4 across the profiles at
# each canonical activation timepoint (0h, LA, 16h, 40h, 5d).
coloc_core <- coloc_out[gene %in% core9_genes]
coloc_core[, gene_symbol := coalesce(sym_lup[gene], gene)]
coloc_core[, tp := parse_timepoint(profile)]
coloc_core <- coloc_core[tp != "other"]
coloc_heat <- coloc_core[, .(max_PP4 = max(PP.H4, na.rm = TRUE)),
                         by = .(gene_symbol, tp)]
coloc_heat[, gene_symbol := factor(gene_symbol, levels = sort(unique(gene_symbol)))]
coloc_heat[, tp := factor(tp, levels = c("0h", "LA", "16h", "40h", "5d"))]

# =============================================================================
# 3. SMR + HEIDI validation
# =============================================================================
# Summary-data-based Mendelian randomization (Zhu et al. 2016, Nat Genet) uses
# only the top cis-eQTL SNP per gene: beta_SMR = beta_GWAS / beta_eQTL, with a
# delta-method standard error. HEIDI then tests whether a single causal
# variant explains the association across the remaining cis-SNPs (chi-square
# on d_i = beta_GWAS_i - beta_SMR * beta_eQTL_i). Dual criteria:
#   p_SMR < 0.05  AND  p_HEIDI > 0.05
# In the current data 70 of the 129 MR-significant genes pass both.

# Best profile per MR-significant gene (IVW or Wald ratio estimate).
mr_sig  <- mr_full[fdr < 0.05 & method %in% c("Inverse variance weighted", "Wald ratio")]
mr_best <- mr_sig[order(fdr, pval), .SD[1], by = gene]

# GWAS prepared for the chr + pos merge and the allele-alignment checks.
gwas_smr <- copy(gwas_raw)
setnames(gwas_smr,
         old = c("chromosome", "base_pair_location", "effect_allele",
                 "other_allele", "beta", "standard_error"),
         new = c("chr_gwas", "pos_gwas", "ea", "oa", "beta_gwas", "se_gwas"))
gwas_smr <- gwas_smr[chr_gwas %in% unique(eqtl_full[gene_id %in% mr_best$gene, chr])]
gwas_smr[, chr_gwas := as.character(chr_gwas)]
setkey(gwas_smr, chr_gwas, pos_gwas)

# Group genes by profile so that each parquet file is read once.
profile_genes <- mr_best[, .(genes = list(gene)), by = profile]
all_results <- list()
result_idx  <- 0L
for (i in seq_len(nrow(profile_genes))) {
  prof         <- profile_genes$profile[i]
  target_genes <- profile_genes$genes[[i]]
  parquet_file <- file.path(parquet_dir, paste0(prof, "_500kb_combined.parquet"))
  if (!file.exists(parquet_file)) next
  # (elided: per-profile parquet loading and variant_id -> chr/pos/ref/alt
  #  parsing -- same pattern as the colocalization worker, retaining
  #  pval_nominal for top-SNP selection)
  for (g in target_genes) {
    gd <- snp_data[phenotype_id == g]
    if (nrow(gd) < 5) next   # insufficient cis coverage -> skip gene
    # Merge with the GWAS by chr + pos.
    setkey(gd, var_chr, var_pos)
    gd <- gwas_smr[gd, nomatch = 0]
    # (elided: identical low-coverage skip guard -- nrow(gd) < 5 -> next)
    # ---- Allele alignment --------------------------------------------------
    # The eQTL slope is the effect of var_alt; the GWAS beta is the effect of
    # ea. Flip the GWAS beta when REF = ea / ALT = oa, keep only aligned or
    # flipped variants, and drop unresolved ones.
    gd[, aligned := (var_alt == ea & var_ref == oa)]
    gd[, flip    := (var_ref == ea & var_alt == oa)]
    gd[flip == TRUE, beta_gwas := -beta_gwas]
    gd <- gd[aligned == TRUE | flip == TRUE]
    # (elided: identical low-coverage skip guard -- nrow(gd) < 5 -> next)
    # ---- SMR test on the top cis-eQTL SNP ----------------------------------
    top_idx  <- which.min(gd$pval_nominal)
    top      <- gd[top_idx]
    beta_eqtl <- top$slope
    se_eqtl   <- top$slope_se
    beta_gwas <- top$beta_gwas
    se_gwas   <- top$se_gwas
    if (abs(beta_eqtl) < 1e-6) next   # near-zero eQTL effect -> skip
    beta_smr <- beta_gwas / beta_eqtl
    se_smr   <- sqrt(se_gwas^2 / beta_eqtl^2 +
                     beta_gwas^2 * se_eqtl^2 / beta_eqtl^4)
    p_smr    <- 2 * pnorm(-abs(beta_smr / se_smr))
    # ---- HEIDI test on the remaining cis-SNPs ------------------------------
    heidi_snps <- gd[-top_idx][abs(slope) > 1e-6]
    k_heidi    <- nrow(heidi_snps)
    if (k_heidi >= 3) {
      # Variance without the covariance term (standard SMR approximation).
      d_i     <- heidi_snps$beta_gwas - beta_smr * heidi_snps$slope
      var_d_i <- heidi_snps$se_gwas^2 + beta_smr^2 * heidi_snps$slope_se^2
      chi_sq  <- sum(d_i^2 / var_d_i)
      p_heidi <- pchisq(chi_sq, k_heidi - 1, lower.tail = FALSE)
    } else {
      p_heidi <- NA_real_   # fewer than 3 HEIDI SNPs -> not evaluable
    }
    result_idx <- result_idx + 1L
    all_results[[result_idx]] <- data.table(
      gene = g, profile = prof, top_snp_rsid = top$rsid,
      top_snp_var_id = top$variant_id,
      beta_eqtl = beta_eqtl, se_eqtl = se_eqtl, pval_eqtl = top$pval_nominal,
      beta_gwas = beta_gwas, se_gwas = se_gwas,
      beta_smr = beta_smr, se_smr = se_smr, p_smr = p_smr,
      n_cis_snps = nrow(gd), n_heidi_snps = k_heidi,
      chi2_heidi = chi_sq, p_heidi = p_heidi)
  }
}

# Combine and flag dual-pass genes (Supplementary Table S5 in the paper).
smr_out <- rbindlist(all_results, fill = TRUE)
smr_out[, pass_smr   := p_smr < 0.05]
smr_out[, pass_heidi := p_heidi > 0.05]
smr_out[is.na(p_heidi), pass_heidi := NA]
smr_out[, pass_both  := pass_smr & pass_heidi]
smr_out[is.na(pass_heidi), pass_both := NA]
setorder(smr_out, p_smr)
n_pass  <- sum(smr_out$pass_both, na.rm = TRUE)   # 70 of 129 in the current data
n_total <- nrow(smr_out)
cat("SMR + HEIDI pass:", n_pass, "of", n_total, "genes\n")
fwrite(smr_out, file.path(RES, "A08_SMR_HEIDI.tsv"), sep = "\t")

# =============================================================================
# 4. Paper Figure 3
# =============================================================================
# (a) Top 25 genes ranked by colocalization support
# (b) Colocalization heatmap: 9 core genes x 5 activation timepoints
# (c) MR vs SMR scatter with SMR/HEIDI validation status

# ---- Figure 3a: top 25 coloc genes bar chart ---------------------------------
coloc_gene_max <- coloc_out[, .(
  max_ratio = max(coloc_ratio, na.rm = TRUE),
  max_PPH4  = max(PP.H4, na.rm = TRUE)), by = gene]
coloc_gene_max <- coloc_gene_max[order(-max_PPH4)][1:25]
coloc_gene_max[, gene_symbol := coalesce(sym_lup[gene], gene)]
coloc_gene_max[, is_core := gene %in% core9_genes]
coloc_gene_max[, gene_symbol := fct_reorder(gene_symbol, max_PPH4)]
fig3a <- ggplot(coloc_gene_max, aes(x = max_PPH4, y = gene_symbol, fill = is_core)) +
  geom_col(width = 0.75) +
  geom_text(aes(x = max_PPH4 + 0.02, label = sprintf("%.3f", max_PPH4)),
            size = 2.5, color = "grey40", hjust = 0) +
  scale_fill_manual(values = c("TRUE" = col_core, "FALSE" = "#AEC6E8"),
                    labels = c("TRUE" = "Core (3 databases)", "FALSE" = "CD4 eQTL only"),
                    name = "") +
  scale_x_continuous(limits = c(0, 1.15), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Top 25 genes with strongest colocalization evidence",
       subtitle = sprintf("%d of %d gene-profile pairs coloc-significant (ratio > 0.7, %d unique genes)",
                          n_coloc, nrow(coloc_out), uniqueN(coloc_out[coloc_sig == TRUE, gene])),
       x = expression("Maximum PP.H"[4] * " (posterior probability of shared causal variant)"),
       y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "bottom", panel.grid.major.y = element_blank())
ggsave(file.path(FIG, "Fig3a_coloc_top25.pdf"), fig3a, width = 6, height = 6.5,
       device = "pdf")

# ---- Figure 3b: coloc heatmap, 9 core genes x 5 timepoints -------------------
fig3b <- ggplot(coloc_heat, aes(x = tp, y = gene_symbol, fill = max_PP4)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(max_PP4 > 0.3, sprintf("%.2f", max_PP4), "")),
            size = 2.5, color = "white") +
  scale_fill_gradientn(colors = c("grey92", "#FFF3CD", "#FF7F0E", "#D62728"),
                       values = c(0, 0.3, 0.7, 1), limits = c(0, 1),
                       name = "PP.H4\n(max)") +
  labs(title = "Colocalization of core AR genes across CD4+ T cell activation",
       subtitle = "Maximum PP.H4 per timepoint | coloc significance: PP.H4/(PP.H3+PP.H4) > 0.7",
       x = "Activation timepoint", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major = element_blank())
ggsave(file.path(FIG, "Fig3b_coloc_heatmap.pdf"), fig3b, width = 6.5, height = 4.5,
       device = "pdf")

# ---- Figure 3c: MR vs SMR scatter with SMR/HEIDI status ----------------------
smr_plot <- merge(smr_out, mr_best[, .(gene, mr_beta = beta)],
                  by = "gene", all.x = TRUE)
smr_plot[, gene_symbol := coalesce(sym_lup[gene], gene)]
smr_plot[, smr_status := fcase(
  pass_smr & pass_heidi,    "SMR + HEIDI pass",
  pass_smr & !pass_heidi,   "SMR pass only",
  default = "SMR not pass")]
smr_plot[, smr_status := factor(smr_status,
  levels = c("SMR + HEIDI pass", "SMR pass only", "SMR not pass"))]
# Label the top 20 SMR-pass genes, split by beta sign for outward placement.
smr_top <- smr_plot[pass_smr == TRUE][order(p_smr)][1:20]
smr_top[, label_side := ifelse(mr_beta > 0, "right", "left")]
fig3c <- ggplot(smr_plot, aes(x = mr_beta, y = beta_smr)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50",
              linewidth = 0.5) +
  geom_point(aes(color = smr_status, size = n_heidi_snps), alpha = 0.75) +
  geom_text_repel(data = smr_top[label_side == "right"],
                  aes(label = gene_symbol),
                  size = 2.5, fontface = "italic",
                  nudge_x = 0.08, direction = "y", hjust = 0,
                  segment.size = 0.3, segment.color = "grey50",
                  box.padding = 0.4, point.padding = 0.3,
                  max.overlaps = 15, min.segment.length = 0.1) +
  # (elided: mirrored geom_text_repel layer for negative-beta genes --
  #  identical to the right-nudge layer above with hjust and nudge_x reversed)
  scale_color_manual(values = c("SMR + HEIDI pass" = col_sig,
                                "SMR pass only" = col_risk,
                                "SMR not pass" = "grey80"),
                     name = "Validation status") +
  scale_size_continuous(range = c(1.5, 4.5), name = "HEIDI SNPs") +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = 0.10)) +
  labs(title = "SMR + HEIDI validation of MR-significant genes",
       subtitle = sprintf("%d/%d genes (%.1f%%) pass SMR+HEIDI | %d SMR-only",
                          n_pass, n_total, 100 * n_pass / n_total,
                          sum(smr_plot$pass_smr & !smr_plot$pass_heidi, na.rm = TRUE)),
       x = expression("MR effect size (" * beta * ", from IVW / Wald ratio)"),
       y = expression("SMR effect size (" * beta[SMR] * ")")) +
  theme_nature() +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "Fig3c_smr_heidi.pdf"), fig3c, width = 8, height = 7,
       device = "pdf")
