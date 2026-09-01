# =============================================================================
# 01_primary_mr.R
# Results 2.1: Transcriptome-wide MR identified a broad allergic rhinitis
# gene landscape in activated CD4+ T cells (paper Figures 1-2)
# Paper: Genetic prioritization of CD4+ T cell-expressed candidate genes for
#        allergic rhinitis (Biomedicines)
# Merged and sanitized from: MR_01_extract_top_eqtl.R, MR_02_map_rsid.R,
#        MR_03_harmonize.R, MR_04_clump.R, MR_05_mr_analysis.R,
#        MR_figures_v4.R, _fig1_workflow.R
# =============================================================================
# Notes:
#  - Run from the repository root; each script assumes cwd == repo root.
#  - Raw GWAS/eQTL data are obtained from public repositories (see the
#    paper's Data Availability statement); they are not distributed here.
#  - Repetitive blocks are elided and marked "# (elided: ...)".
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries, paths and shared constants
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)     # fork-based parallelism (mclapply)
  library(TwoSampleMR)  # harmonisation, clumping, Steiger, MR estimation
  library(mr.raps)      # robust adjusted profile score
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
})
# Repository-relative paths. Raw public data (see Data Availability) are
# expected under data/ (Soskic_CD4/, snp150_by_chr/, gwas/, g1000_eur/, plink/).
DATA <- "data"
RES  <- "results"
FIG  <- "figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
plink_bin <- normalizePath(file.path(DATA, "plink/plink"))
bfile     <- normalizePath(file.path(DATA, "g1000_eur/g1000_eur"))
# Five canonical activation time points (Soskic et al., 2022)
tp_levels <- c("0h", "LA", "16h", "40h", "5d")
# Nine core genes replicated in three eQTL databases (Soskic, DICE, eQTLGen)
core9_symbols <- c("NFKB1", "IL18RAP", "AHI1", "IKZF3",
                   "SLC25A46", "TEF", "TLR1", "ZBTB38", "CD247")

# -----------------------------------------------------------------------------
# 1. Extract top cis-eQTL per gene-profile pair (Soskic CD4+ T cells)
# -----------------------------------------------------------------------------
# One tensorQTL cis-eQTL file (variants within a 500 kb window per gene) per profile.
tensor_dir <- file.path(DATA, "Soskic_CD4")
files <- list.files(tensor_dir, pattern = "_500kb_window_tensorQTL\\.txt$",
                    full.names = TRUE)
# (elided: 46 individual per-profile file reads; the generic loop below
#  applies the identical read + profile-tagging logic to every file)
all_profiles <- lapply(files, function(f) {
  dt <- fread(f)
  dt[, profile := sub("_500kb_window_tensorQTL\\.txt$", "", basename(f))]
  dt
})
eqtl <- rbindlist(all_profiles)
# Instrument selection: permutation FDR < 0.05, F >= 10 (F = (slope/slope_se)^2)
eqtl <- eqtl[pval_beta < 0.05]
eqtl[, F_stat := (slope / slope_se)^2]
eqtl <- eqtl[F_stat >= 10]
eqtl_out <- eqtl[, .(
  profile,
  gene_id       = phenotype_id,
  snp_id        = variant_id,
  chr           = sub("_.*", "", variant_id),
  pos           = as.integer(sub("^[^_]+_([^_]+)_.*", "\\1", variant_id)),
  ref           = sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", variant_id),
  alt           = sub("^.*_([^_]+)$", "\\1", variant_id),
  beta          = slope,
  se            = slope_se,
  pval_nominal,
  pval_beta,
  maf,
  F_stat
)]
fwrite(eqtl_out, file.path(RES, "01_top_eqtl.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 2. Map variant positions to rsIDs (UCSC snp150, per-chromosome files)
# -----------------------------------------------------------------------------
# The full snp150 table is ~7.1 GB; per-chromosome sub-files are used instead.
eqtl <- fread(file.path(RES, "01_top_eqtl.tsv"))
eqtl[, chrom := paste0("chr", chr)]
pos_by_chr <- split(eqtl$pos, eqtl$chrom)
chrs       <- names(pos_by_chr)
split_dir  <- file.path(DATA, "snp150_by_chr")
# (elided: 22 per-chromosome snp150 lookups; the mclapply call below applies
#  the same position -> rsID filter to each chromosome file)
snp_ref <- rbindlist(mclapply(chrs, function(ch) {
  f <- file.path(split_dir, paste0(ch, ".txt"))
  if (!file.exists(f)) return(NULL)
  pos_needed <- pos_by_chr[[ch]]
  dt <- fread(f, select = c(2, 4, 5), col.names = c("chrom", "pos", "rsid"))
  dt[chrom == ch & pos %in% pos_needed]
}, mc.cores = min(12, length(chrs))))
snp_ref[, chr := sub("chr", "", chrom)]
eqtl[, chr := as.character(chr)]
eqtl <- merge(eqtl, snp_ref[, .(chr, pos, rsid)], by = c("chr", "pos"), all.x = TRUE)
eqtl <- eqtl[!is.na(rsid)]
eqtl[, chrom := NULL]
fwrite(eqtl, file.path(RES, "02_eqtl_with_rsid.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 3. MHC exclusion and harmonisation against the AR GWAS
# -----------------------------------------------------------------------------
# MHC region chr6:25.5-34 Mb excluded (complex long-range LD violates MR assumptions).
eqtl <- eqtl[!(chr == 6 & pos >= 25500000 & pos <= 34000000)]
# GCST90468131: hayfever, UK Biobank Quickdraws, GRCh38, N = 394,626, case
# fraction ~20%. Standard GWAS loading pattern used throughout the pipeline.
# Note: a data.table J() subset drops the key -- call setkey() before further joins.
gwas_raw <- fread(
  cmd = paste("gunzip -c", shQuote(file.path(DATA, "gwas/GCST90468131.h.tsv.gz"))),
  sep = "\t"
)
gwas_raw <- gwas_raw[!is.na(beta) & !is.na(standard_error) & standard_error > 0]
gwas_raw <- gwas_raw[!is.na(rsid) & rsid != ""]
setkey(gwas_raw, rsid)
gwas <- gwas_raw[, .(
  SNP                   = rsid,
  beta.outcome          = beta,
  se.outcome            = standard_error,
  pval.outcome          = p_value,
  eaf.outcome           = effect_allele_frequency,
  effect_allele.outcome = effect_allele,
  other_allele.outcome  = other_allele,
  outcome               = "GCST90468131",
  id.outcome            = "GCST90468131"
)]
gwas <- gwas[SNP %in% eqtl$rsid]
exposure <- eqtl[, .(
  SNP                    = rsid,
  beta.exposure          = beta,
  se.exposure            = se,
  pval.exposure          = pval_nominal,
  eaf.exposure           = maf,
  effect_allele.exposure = alt,
  other_allele.exposure  = ref,
  exposure               = gene_id,
  id.exposure            = gene_id,
  profile                = profile,
  F_stat                 = F_stat
)]
# action = 2: effect-allele alignment; palindromic SNPs dropped if strand ambiguous.
dat <- harmonise_data(
  exposure_dat = as.data.frame(exposure),
  outcome_dat  = as.data.frame(gwas),
  action       = 2
)
dat <- dat[dat$mr_keep == TRUE, ]
fwrite(dat, file.path(RES, "03_harmonized.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 4. LD clumping and Steiger directionality filtering
# -----------------------------------------------------------------------------
dat <- fread(file.path(RES, "03_harmonized.tsv"))
profiles <- unique(dat$profile)
# Clumping per profile (PLINK v1.9, r2 < 0.001, 500 kb, 1000 Genomes EUR).
# (elided: repeated per-profile clump calls; one representative instance of
#  the per-profile clump + SNP retention logic is shown below)
clumped_list <- lapply(profiles, function(p) {
  sub_dat <- dat[profile == p]
  exposure_for_clump <- data.frame(
    SNP           = sub_dat$SNP,
    beta.exposure = sub_dat$beta.exposure,
    se.exposure   = sub_dat$se.exposure,
    pval.exposure = sub_dat$pval.exposure,
    eaf.exposure  = sub_dat$eaf.exposure,
    effect_allele.exposure = sub_dat$effect_allele.exposure,
    other_allele.exposure  = sub_dat$other_allele.exposure,
    exposure      = sub_dat$exposure,
    id.exposure   = sub_dat$exposure
  )
  clumped <- tryCatch(
    clump_data(
      exposure_for_clump,
      clump_r2   = 0.001,
      clump_kb   = 500,
      pop        = "EUR",
      bfile      = bfile,
      plink_bin  = plink_bin
    ),
    error = function(e) NULL
  )
  if (is.null(clumped) || nrow(clumped) == 0) return(NULL)
  sub_dat[SNP %in% clumped$SNP]
})
dat_clumped <- rbindlist(clumped_list, fill = TRUE)
dat_clumped <- dat_clumped[F_stat >= 10]  # re-apply F threshold after clumping
# Steiger test: variant explains more exposure variance than outcome (N_exp = 119, N_out = 394,626).
dat_clumped$samplesize.exposure <- rep(119,    nrow(dat_clumped))
dat_clumped$samplesize.outcome  <- rep(394626, nrow(dat_clumped))
dat_clumped <- steiger_filtering(as.data.frame(dat_clumped))
dat_clumped <- dat_clumped[dat_clumped$steiger_dir == TRUE |
                           is.na(dat_clumped$steiger_dir), ]
fwrite(dat_clumped, file.path(RES, "04_clumped.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 5. Two-sample MR per gene-profile pair (Wald ratio / IVW / weighted median /
#    weighted mode / MR-Egger / MR-RAPS) with BH FDR correction
# -----------------------------------------------------------------------------
dat <- as.data.frame(fread(file.path(RES, "04_clumped.tsv")))
genes    <- unique(dat$exposure)
profiles <- unique(dat$profile)
results_list <- list()
for (gene in genes) {
  for (p in profiles) {
    sub <- dat[dat$exposure == gene & dat$profile == p, ]
    if (nrow(sub) == 0) next
    # Single instrument after clumping: Wald ratio estimator
    if (nrow(sub) == 1) {
      wr <- mr(sub, method_list = "mr_wald_ratio")
      results_list[[length(results_list) + 1]] <- data.frame(
        gene = gene, profile = p, nsnp = 1, method = "Wald ratio",
        beta = wr$b, se = wr$se, pval = wr$pval)
      next
    }
    # Multiple instruments: IVW plus pleiotropy-robust sensitivity analyses
    mr_res <- mr(sub, method_list = c(
      "mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"
    ))
    raps_res <- tryCatch(
      mr.raps(
        b_exp  = sub$beta.exposure,
        b_out  = sub$beta.outcome,
        se_exp = sub$se.exposure,
        se_out = sub$se.outcome,
        over.dispersion = TRUE,   # systematic + idiosyncratic pleiotropy
        loss.function   = "huber" # robustness to outlier instruments
      ),
      error = function(e) NULL
    )
    rows <- data.frame(
      gene = gene, profile = p, nsnp = mr_res$nsnp, method = mr_res$method,
      beta = mr_res$b, se = mr_res$se, pval = mr_res$pval)
    if (!is.null(raps_res)) {
      rows <- rbind(rows, data.frame(
        gene = gene, profile = p, nsnp = nrow(sub), method = "MR-RAPS",
        beta = raps_res$beta.hat, se = raps_res$beta.se,
        pval = 2 * pnorm(-abs(raps_res$beta.hat / raps_res$beta.se))))
    }
    results_list[[length(results_list) + 1]] <- rows
  }
}
results <- rbindlist(results_list, fill = TRUE)
# BH FDR over all tested pairs; IVW and Wald ratio form the primary results.
ivw_results <- results[results$method %in% c("Inverse variance weighted",
                                             "Wald ratio"), ]
ivw_results$fdr <- p.adjust(ivw_results$pval, method = "BH")
sig_genes <- unique(ivw_results$gene[ivw_results$fdr < 0.05])
fwrite(results,     file.path(RES, "05_MR_all_methods.tsv"), sep = "\t")
fwrite(ivw_results, file.path(RES, "05_MR_ivw_fdr.tsv"),    sep = "\t")
# Expected counts: 50,654 pairs across 10,328 genes scanned; 129 genes at FDR < 0.05.

# -----------------------------------------------------------------------------
# 6. Shared figure data, theme and helpers
# -----------------------------------------------------------------------------
# Gene symbol mapping (Ensembl ID -> HGNC symbol) from an annotation helper script.
sym_dt <- fread(file.path(RES, "gene_symbol_full.tsv"))
setnames(sym_dt, c("gene", "gene_symbol"))
sym_lup <- setNames(sym_dt$gene_symbol, sym_dt$gene)
# Coloc results (Results 2.2); used for the colocalization aesthetic of Fig 2c.
coloc <- fread(file.path(RES, "06_coloc.tsv"))
# Map the 46 expression profiles to the five canonical activation time points
parse_timepoint <- function(profile) {
  tp <- case_when(
    grepl("(?<!4)0h", profile, perl = TRUE) ~ "0h",
    grepl("_LA$|_LA_", profile) ~ "LA",
    grepl("16h",  profile) ~ "16h",
    grepl("40h",  profile) ~ "40h",
    grepl("5d",   profile) ~ "5d",
    TRUE ~ "other"
  )
  factor(tp, levels = c(tp_levels, "other"))
}
# Best-associated profile per significant gene (lowest p-value; FDR < 0.05)
mr_full <- fread(file.path(RES, "05_MR_ivw_fdr.tsv"))
mr_sig  <- mr_full[fdr < 0.05]
mr_best <- mr_sig[order(fdr, pval)][!duplicated(gene)]
mr_best[, `:=`(
  timepoint   = parse_timepoint(profile),
  direction   = ifelse(beta > 0, "Risk-increasing", "Protective"),
  gene_symbol = coalesce(sym_lup[gene], gene),
  log10fdr    = -log10(fdr)
)]
mr_best[, is_core := gene_symbol %in% core9_symbols]
n_sig  <- nrow(mr_best)                 # 129 genes
n_core <- length(core9_symbols)         # 9 core genes
# Shared plotting theme (high-contrast, visible grid lines, PDF only)
theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
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
      plot.margin      = margin(6, 6, 6, 6)
    )
}
# High-contrast NPG palette
col_risk    <- "#D62728"   # risk-increasing
col_protect <- "#1F77B4"   # protective

# -----------------------------------------------------------------------------
# 7. Figure 1 - study design workflow (paper Figure 1)
# -----------------------------------------------------------------------------
# Manual schematic of the three-module design; Module A is drawn to the point
# covered by this section. Downstream coloc / SMR+HEIDI / timepoint steps and
# per-step boxes of Modules B and C are elided (one endpoint box per module).
# (elided: downstream Module A boxes and per-step Module B/C boxes -- one
#  representative box per module is kept below)
make_box   <- function(x, y, w, h, label, fill, text_col = "#FFFFFF") list(
  x = x, y = y, w = w, h = h, label = label, fill = fill, text_col = text_col)
make_arrow <- function(x0, y0, x1, y1) list(
  x0 = x0, y0 = y0, x1 = x1, y1 = y1)
x_start <- 1
y_a <- 7.5   # Module A (primary MR)
y_b <- 5.0   # Module B (replication)
y_c <- 2.5   # Module C (mediation)
box_h <- 1.4
a_boxes <- list(
  make_box(x_start,        y_a, 2.2, box_h, "Soskic CD4+ eQTL\n46 profiles · 5 timepoints", "#4C72B0"),
  make_box(x_start +  2.8, y_a, 2.2, box_h, "Instrument Selection\nF >= 10 · LD clump · Steiger", "#4C72B0"),
  make_box(x_start +  5.6, y_a, 1.8, box_h, "MHC Exclusion\nchr6:25.5-34Mb", "#7B9DC8"),
  make_box(x_start +  8.0, y_a, 2.0, box_h, "MR Analysis\nIVW · Wald · RAPS", "#4C72B0")
)
# (elided: DICE per-cell-type and eQTLGen step boxes of Module B)
b_boxes <- list(
  make_box(x_start,        y_b, 2.2, box_h, "DICE + eQTLGen eQTL\n15 cell types + whole blood", "#5A9E6F"),
  make_box(x_start + 11.0, y_b, 2.2, box_h, "Cross-DB Comparison\nVenn · 9 core genes", "#70B080")
)
# (elided: two-step metabolite MR and coloc boxes of Module C)
c_boxes <- list(
  make_box(x_start, y_c, 2.5, box_h, "Gene -> Metabolite MR\n150 genes -> 243 metabolites", "#D4954B")
)
gwas_box <- list(x = x_start + 4.2, y = 8.8, w = 2.0, h = 0.9,
                 label = "AR GWAS\nGCST90468131 · N=394,626",
                 fill = "#E8E0D8", text_col = "#444444")
box_dt <- rbindlist(lapply(c(a_boxes, b_boxes, c_boxes), as.data.table))
box_dt <- rbind(box_dt, as.data.table(gwas_box), fill = TRUE)
a_arrows <- list(
  make_arrow(3.2, y_a, 2.8, y_a), make_arrow(5.0, y_a, 5.6, y_a),
  make_arrow(7.4, y_a, 8.0, y_a)
)
# (elided: remaining Module A arrows connecting the downstream boxes)
b_arrows <- list(make_arrow(3.2, y_b, 3.0, y_b),
                 make_arrow(5.0, y_b, 5.6, y_b),
                 make_arrow(10.2, y_b, 11.0, y_b))
arrow_dt <- rbindlist(lapply(c(a_arrows, b_arrows), as.data.table))
gwas_arrow <- data.table(x0 = 5.2, y0 = 8.8, x1 = 5.2, y1 = 8.2)
mod_labels <- data.table(
  x = 0.3, y = c(y_a, y_b, y_c),
  label = c("Module A\nPrimary MR", "Module B\nCross-DB Replication",
            "Module C\nMetabolite Mediation"),
  color = c("#4C72B0", "#5A9E6F", "#D4954B")
)
legend_dt <- data.table(
  x = c(0.8, 7, 13), y = 0.6,
  label = c("AR GWAS outcome (all modules share GCST90468131, N = 394,626, case fraction ~20%)",
            "MHC region excluded: chr6:25.5-34.0 Mb  |  Soskic N = 119 European donors",
            "5 activation timepoints: 0h / LA / 16h / 40h / 5d"),
  color = c("#666666", "#888888", "#888888")
)
p_fig1 <- ggplot() +
  geom_rect(data = box_dt,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
                fill = fill), color = NA, alpha = 0.92) +
  geom_text(data = box_dt, aes(x = x, y = y, label = label, color = text_col),
            size = 2.5, lineheight = 0.9, fontface = "bold") +
  geom_segment(data = arrow_dt, aes(x = x0 + 0.08, y = y0, xend = x1 - 0.08, yend = y1),
               arrow = arrow(length = unit(2, "pt"), type = "closed"),
               linewidth = 0.5, color = "#888888") +
  geom_segment(data = gwas_arrow, aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(2, "pt"), type = "closed"),
               linewidth = 0.5, color = "#999999", linetype = "dotted") +
  geom_text(data = mod_labels, aes(x = x, y = y, label = label, color = color),
            size = 3.2, fontface = "bold", hjust = 0, lineheight = 0.9) +
  geom_text(data = legend_dt, aes(x = x, y = y, label = label, color = color),
            size = 2.2, hjust = 0) +
  scale_color_identity() +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 18.5), ylim = c(0, 9.2)) +
  theme_void()
# Export as PDF (185 x 90 mm)
pdf(file.path(FIG, "Fig1_workflow.pdf"), width = 185 / 25.4, height = 90 / 25.4,
    useDingbats = FALSE)
print(p_fig1)
dev.off()

# -----------------------------------------------------------------------------
# 8. Figure 2a - genomic landscape of MR-significant genes (paper Figure 2a)
# -----------------------------------------------------------------------------
# Manhattan-style layout of the best per-gene MR result across autosomes after
# MHC exclusion; positions come from the top eQTL variant.
# (Panel code adapted from the original figure scripts.)
# Best profile per gene across ALL genes (background = non-significant)
man_best <- mr_full[order(pval), .SD[1], by = gene]
man_best[, `:=`(gene_symbol = coalesce(sym_lup[gene], gene),
                log10fdr    = -log10(pmin(fdr, 1)),
                is_core     = gene_symbol %in% core9_symbols)]
gene_pos <- fread(file.path(RES, "01_top_eqtl.tsv"),
                  select = c("gene_id", "chr", "pos"))
gene_pos <- unique(gene_pos, by = "gene_id")[
  , .(chr = as.integer(chr[1]), gene_pos = as.integer(median(pos))),
  by = .(gene = gene_id)]
man_plot <- merge(man_best, gene_pos, by = "gene", all.x = TRUE)
man_plot <- man_plot[!is.na(chr) & chr <= 22][order(chr, gene_pos)]
chr_sizes <- man_plot[, .(max_pos = max(gene_pos, na.rm = TRUE)), by = chr][order(chr)]
chr_sizes[, offset := c(0, cumsum(as.numeric(max_pos))[1:(.N - 1)]) + 5e7 * (chr - 1)]
man_plot <- merge(man_plot, chr_sizes[, .(chr, offset)], by = "chr", all.x = TRUE)
man_plot[, x_pos := gene_pos + offset]
chr_lab <- man_plot[, .(mid_x = mean(x_pos, na.rm = TRUE)), by = chr][order(chr)]
# Alternating chromosome shades; significant genes colored, the rest grey
man_plot[, `:=`(is_sig   = fdr < 0.05,
                fill_col = fifelse(chr %% 2 == 0, "#7E9BBF", "#3C5488"))]
# MHC region band (chr6:25.5-34 Mb, coordinates mapped to the cumulative axis)
chr6_offset <- chr_sizes[chr == 6, offset]
mhc_start <- chr6_offset + 25.5e6
mhc_end   <- chr6_offset + 34.0e6
# Labels: top 25 significant genes by FDR plus the 9 core genes
top25_man <- man_plot[is_sig == TRUE][order(fdr)][1:25, gene]
man_plot[, label_m := ifelse(gene %in% union(top25_man, man_best[is_core == TRUE, gene]) &
                               is_sig == TRUE,
                             gene_symbol, NA_character_)]
p_fig2a <- ggplot(man_plot, aes(x = x_pos, y = log10fdr)) +
  annotate("rect", xmin = mhc_start, xmax = mhc_end, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +
  annotate("text", x = (mhc_start + mhc_end) / 2, y = max(man_plot$log10fdr) * 1.03,
           label = "MHC\nexcluded", size = 2.4, color = "grey50", lineheight = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  geom_point(data = man_plot[is_sig == FALSE],
             color = "grey75", size = 0.7, alpha = 0.5, shape = 16) +
  geom_point(data = man_plot[is_sig == TRUE & is_core == FALSE],
             aes(color = fill_col), size = 2.0, alpha = 0.85, shape = 16) +
  geom_point(data = man_plot[is_core == TRUE],
             aes(color = fill_col), size = 3.5, alpha = 1, shape = 18, stroke = 0.3) +
  scale_color_identity() +
  geom_text_repel(aes(label = label_m),
                  size = 2.3, fontface = "italic", na.rm = TRUE,
                  max.overlaps = 35, segment.size = 0.3, segment.color = "grey50",
                  box.padding = 0.3, point.padding = 0.15,
                  direction = "y", nudge_y = 0.2) +
  geom_vline(xintercept = chr_sizes$offset[-1] - 2.5e7,
             color = "grey85", linewidth = 0.3) +
  scale_x_continuous(breaks = chr_lab$mid_x, labels = chr_lab$chr,
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title = "Genomic landscape of MR-significant genes across the autosomes",
    subtitle = sprintf("Best profile per gene (n = %d significant at FDR < 0.05) | diamonds = %d core genes | MHC region (chr6:25.5-34 Mb) excluded",
                       n_sig, n_core),
    x = "Chromosome (GRCh38)",
    y = expression(-log[10]("FDR"))
  ) +
  theme_nature() +
  theme(panel.grid.major = element_blank())
ggsave(file.path(FIG, "Fig2a_landscape.pdf"), p_fig2a, width = 14, height = 5,
       device = "pdf")

# -----------------------------------------------------------------------------
# 9. Figure 2b - volcano plot of the best MR profile per gene (paper Figure 2b)
# -----------------------------------------------------------------------------
mr_best_plot <- copy(mr_best)
# Labels: 9 core genes + top 15 genes by FDR
label_set <- union(core9_symbols, mr_best_plot[order(fdr)][1:15, gene_symbol])
mr_best_plot[, label_gene := ifelse(gene_symbol %in% label_set,
                                    gene_symbol, NA_character_)]
p_fig2b <- ggplot(mr_best_plot, aes(x = beta, y = log10fdr)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50", linewidth = 0.5) +
  geom_point(data = mr_best_plot[is_core == FALSE],
             aes(color = direction, size = log10fdr), alpha = 0.65, stroke = 0.2) +
  geom_point(data = mr_best_plot[is_core == TRUE],
             aes(color = direction, size = log10fdr + 1.5), alpha = 0.95,
             stroke = 0.5, shape = 18) +
  geom_text_repel(
    aes(label = label_gene),
    size = 2.4, max.overlaps = 25, fontface = "italic",
    segment.size = 0.3, segment.color = "grey50",
    box.padding = 0.3, point.padding = 0.2, na.rm = TRUE
  ) +
  scale_color_manual(
    values = c("Risk-increasing" = col_risk, "Protective" = col_protect),
    name = "Effect direction"
  ) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(
    title = "CD4+ T cell gene expression MR effects on allergic rhinitis",
    subtitle = sprintf("Best profile per gene (n = %d, FDR < 0.05) | diamonds = %d core genes (3 databases)",
                       n_sig, n_core),
    x = expression("MR effect size (" * beta * ", per SD increase in gene expression)"),
    y = expression(-log[10]("FDR"))
  ) +
  theme_nature() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = NA))
ggsave(file.path(FIG, "Fig2b_volcano.pdf"), p_fig2b, width = 7.5, height = 5.5,
       device = "pdf")

# -----------------------------------------------------------------------------
# 10. Figure 2c - forest plot of the 9 core genes (paper Figure 2c)
# -----------------------------------------------------------------------------
core9_mr <- mr_best[is_core == TRUE]
core9_mr[, `:=`(
  ci_lo = beta - 1.96 * se,
  ci_hi = beta + 1.96 * se,
  gene_symbol = factor(gene_symbol, levels = gene_symbol[order(-abs(beta))])
)]
core9_mr[, coloc_sig := gene %in% coloc[coloc_sig == TRUE, gene]]
p_fig2c <- ggplot(core9_mr, aes(x = beta, y = gene_symbol,
                                color = direction, shape = coloc_sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.35, linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("Risk-increasing" = col_risk, "Protective" = col_protect),
    name = "Direction"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 18),
    labels = c("FALSE" = "No", "TRUE" = "Yes (PP.H4/(PP.H3+PP.H4) > 0.7)"),
    name = "Colocalization"
  ) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = "Replicated in 3 databases: Soskic CD4+ & DICE & eQTLGen",
           size = 2.5, color = "grey40", fontface = "italic") +
  labs(
    title = "Core causal genes replicated across three eQTL databases",
    subtitle = "Allergic rhinitis (n = 9 genes in CD4+ T cell & DICE & eQTLGen)",
    x = expression("MR effect size (" * beta * "  " %+-% "  95% CI, CD4+ eQTL data)"),
    y = NULL
  ) +
  theme_nature() +
  theme(axis.text.y = element_text(face = "italic", size = 10))
ggsave(file.path(FIG, "Fig2c_core_forest.pdf"), p_fig2c, width = 6.5, height = 5,
       device = "pdf")
