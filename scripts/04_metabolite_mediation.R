# =============================================================================
# 04_metabolite_mediation.R
# Results 2.4: Exploratory metabolite screening identified indirect-effect
# signals, but failed the shared-variant test (paper Figure 5)
# Paper: Genetic prioritization of CD4+ T cell-expressed candidate genes for
#        allergic rhinitis (Biomedicines)
# Merged and sanitized from: MR_14_gene_metabolite_mr.R, MR_15_metabolite_ar_mr.R,
#        MR_16_mediation.R, MR_17_metabolite_coloc.R, MR_figures_v4.R
# =============================================================================
# Notes:
#  - Run from the repository root; each script assumes cwd == repo root.
#  - Raw GWAS/metabolite data are obtained from public repositories (see the
#    paper's Data Availability statement); they are not distributed here.
#  - Repetitive blocks are elided and marked "# (elided: ...)".
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Libraries, paths and shared settings
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(mr.raps)
  library(parallel)
  library(coloc)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(forcats)
  library(patchwork)
})
# Raw data come from public repositories (see Data Availability); adjust these
# paths to local copies.
met_dir     <- "data/metabolites"                       # 314 metabolite GWAS (GRCh38)
gwas_file   <- "data/allergic_rhinitis_gwas_data/GCST90468131/GCST90468131.h.tsv.gz"
snp150_file <- "data/snp150_hg38.txt.gz"
plink_bin   <- "data/ld_ref/plink"                      # PLINK binary
bfile       <- "data/ld_ref/g1000_eur"                  # 1000 Genomes EUR reference
res_dir <- "results"
fig_dir <- "figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
n_cores <- 4

# -----------------------------------------------------------------------------
# 1. Gene-to-metabolite MR (paper Figure 5a)
#    Exposures: variants of MR-significant CD4+ T cell eQTL genes (FDR < 0.05
#    in the primary gene -> AR analysis) shared with each metabolite GWAS at
#    the same chr:pos; Wald ratio for a single shared variant, IVW + weighted
#    median otherwise. Expected: 1,693 significant pairs across 126 genes and
#    314 metabolite GWAS (FDR < 0.05).
# -----------------------------------------------------------------------------
eqtl <- fread(file.path(res_dir, "02_eqtl_with_rsid.tsv"))
sig_genes <- unique(fread(file.path(res_dir, "05_MR_ivw_fdr.tsv"))[fdr < 0.05, gene])
eqtl_sig <- eqtl[gene_id %in% sig_genes]
eqtl_sig[, chr_pos := paste0(chr, "_", pos)]
# Shell-level awk pre-filter so each ~340 MB metabolite file is reduced to the
# requested eQTL positions before fread sees it (pattern reused in Section 2)
pos_dt <- unique(eqtl_sig[, .(chr = as.character(chr), pos = as.character(pos))])
pos_file <- tempfile(fileext = ".tsv")
fwrite(pos_dt, pos_file, sep = "\t", col.names = FALSE)
filter_script <- tempfile(fileext = ".sh")
writeLines(c("#!/bin/bash",
             "awk -F'\\t' 'NR==FNR{pos[$1\"\\t\"$2]=1;next} FNR==1||(($1\"\\t\"$2) in pos)' \"$1\" <(gunzip -c \"$2\")"),
           filter_script)
system2("chmod", c("+x", filter_script))
metabolite_files <- list.files(met_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
# (elided: enumeration of all 314 metabolite GWAS files; each is processed by the routine below)
genes <- unique(eqtl_sig$gene_id)
# Representative per-file routine: harmonize the gene eQTL variants against one
# metabolite GWAS and run per-gene MR.
process_metabolite <- function(met_file) {
  met_id <- sub("_buildGRCh38\\.tsv\\.gz$", "", basename(met_file))
  met <- tryCatch(fread(cmd = paste(filter_script, pos_file, met_file), header = TRUE),
                  error = function(e) NULL)
  if (is.null(met) || nrow(met) == 0) return(NULL)
  setnames(met, old = names(met),
           new = c("chromosome", "base_pair_location", "effect_allele",
                   "other_allele", "effect_allele_frequency",
                   "beta", "standard_error", "p_value", "variant_id"))
  met <- met[!is.na(beta) & !is.na(standard_error) & standard_error > 0]
  if (nrow(met) == 0) return(NULL)
  met[, chr_pos := paste0(chromosome, "_", base_pair_location)]
  res_list <- list()
  for (g in genes) {
    eqtl_gene <- eqtl_sig[gene_id == g]
    common <- intersect(eqtl_gene$chr_pos, met$chr_pos)
    if (length(common) == 0) next
    exp_sub <- eqtl_gene[chr_pos %in% common][, .SD[which.min(pval_nominal)], by = chr_pos]
    met_sub <- met[chr_pos %in% common][, .SD[which.min(p_value)], by = chr_pos]
    merged <- merge(exp_sub[, .(chr_pos, rsid, beta, se, pval_nominal, maf, alt, ref)],
                    met_sub[, .(chr_pos, beta, standard_error, p_value,
                                effect_allele_frequency, effect_allele, other_allele)],
                    by = "chr_pos")
    if (nrow(merged) == 0) next
    exposure_df <- data.frame(SNP = merged$rsid, beta.exposure = merged$beta.x,
                              se.exposure = merged$se, pval.exposure = merged$pval_nominal,
                              eaf.exposure = merged$maf,
                              effect_allele.exposure = merged$alt,
                              other_allele.exposure = merged$ref,
                              exposure = g, id.exposure = g)
    outcome_df <- data.frame(SNP = merged$rsid, beta.outcome = merged$beta.y,
                             se.outcome = merged$standard_error,
                             pval.outcome = merged$p_value,
                             eaf.outcome = merged$effect_allele_frequency,
                             effect_allele.outcome = merged$effect_allele,
                             other_allele.outcome = merged$other_allele,
                             outcome = met_id, id.outcome = met_id)
    dat <- harmonise_data(exposure_df, outcome_df, action = 2)
    dat <- dat[dat$mr_keep == TRUE, ]
    if (nrow(dat) == 0) next
    if (nrow(dat) == 1) {
      wr <- mr(dat, method_list = "mr_wald_ratio")
      res_list[[length(res_list) + 1]] <- data.frame(gene = g, metabolite = met_id,
        nsnp = 1, method = "Wald ratio", beta = wr$b, se = wr$se, pval = wr$pval)
    } else {
      mr_res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median"))
      res_list[[length(res_list) + 1]] <- data.frame(gene = g, metabolite = met_id,
        nsnp = mr_res$nsnp, method = mr_res$method, beta = mr_res$b,
        se = mr_res$se, pval = mr_res$pval)
    }
  }
  if (length(res_list) == 0) return(NULL)
  rbindlist(res_list, fill = TRUE)
}
results_list <- mclapply(metabolite_files, process_metabolite, mc.cores = n_cores)
results <- rbindlist(Filter(Negate(is.null), results_list), fill = TRUE)
ivw <- results[results$method %in% c("Inverse variance weighted", "Wald ratio"), ]
ivw$fdr <- p.adjust(ivw$pval, method = "BH")
# Expected: 1,693 significant gene-metabolite pairs (FDR < 0.05)
fwrite(results, file.path(res_dir, "C01_gene_metabolite_MR_all.tsv"), sep = "\t")
fwrite(ivw,     file.path(res_dir, "C01_gene_metabolite_MR_ivw.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 2. Metabolite-to-AR MR (paper Figure 5b)
#    Instruments: genome-wide-significant (p < 5e-8) variants of each
#    C01-significant metabolite, clumped at r2 < 0.001 / 500 kb (EUR). The AR
#    GWAS (~12.7M rows) is never fully loaded: a rsid whitelist is written to a
#    temp file and awk pre-filters the GWAS at the shell level, so fread only
#    sees the whitelisted rows (deferred loading optimization).
# -----------------------------------------------------------------------------
gene_met <- fread(file.path(res_dir, "C01_gene_metabolite_MR_ivw.tsv"),
                  select = c("metabolite", "fdr"))
sig_metabolites <- unique(gene_met[fdr < 0.05]$metabolite)
if (length(sig_metabolites) == 0) stop("No significant metabolites to test")
met_files <- file.path(met_dir, paste0(sig_metabolites, "_buildGRCh38.tsv.gz"))
met_files <- met_files[file.exists(met_files)]
# Phase 1 -- GWS scan: awk keeps only rows with p < 5e-8 per file, so fread
# never sees the full metabolite file.
scan_gws <- function(mf) {
  mid <- sub("_buildGRCh38\\.tsv\\.gz$", "", basename(mf))
  dt <- tryCatch(fread(cmd = sprintf("gunzip -c '%s' | awk -F'\\t' 'NR==1||($8!=\"\"&&$8+0<5e-8)'", mf),
                       select = c("chromosome", "base_pair_location", "effect_allele",
                                  "other_allele", "effect_allele_frequency", "beta",
                                  "standard_error", "p_value")),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt <- dt[!is.na(beta) & !is.na(standard_error) & standard_error > 0]
  dt <- dt[(beta / standard_error)^2 >= 10]
  if (nrow(dt) == 0) return(NULL)
  dt[, met_id := mid]
  dt
}
all_gws <- rbindlist(Filter(Negate(is.null),
                            mclapply(met_files, scan_gws, mc.cores = n_cores)), fill = TRUE)
# (elided: GWS-scan invocations over the remaining significant metabolite files; same mclapply pattern as Section 1)
if (nrow(all_gws) == 0) stop("No genome-wide-significant hits")
# Phase 2 -- resolve chr:pos to rsid via a single-pass snp150 lookup.
# (elided: snp150 chr:pos -> rsid lookup; same temp-file + awk pre-filter
#  pattern as the position filter in Section 1, extracting chr/pos/rsid)
all_gws[, chr := as.character(chromosome)]
all_gws[, pos := as.integer(base_pair_location)]
all_gws <- merge(all_gws, snp_lkp, by = c("chr", "pos"), all.x = TRUE)
all_gws <- all_gws[!is.na(rsid)]
# Phase 2.5 -- AR GWAS, deferred loading: rsid whitelist written to a temp file,
# awk pre-filters the GWAS at the shell level so fread sees only the
# whitelisted rows instead of the full ~12.7M-row file.
gw_rsids <- unique(all_gws$rsid)
rsid_file <- tempfile(fileext = ".tsv")
fwrite(data.table(rsid = gw_rsids), rsid_file, sep = "\t", col.names = FALSE)
gw_filter_sh <- tempfile(fileext = ".sh")
writeLines(c("#!/bin/bash",
             paste0("awk -F'\\t' 'NR==FNR{k[$1]=1;next}",
                    " FNR==1 || ($10 in k) {print $10\"\\t\"$5\"\\t\"$6\"\\t\"$8\"\\t\"$7\"\\t\"$3\"\\t\"$4}'",
                    " \"$1\" <(gunzip -c \"$2\")")),
           gw_filter_sh)
system2("chmod", c("+x", gw_filter_sh))
gwas_raw <- fread(cmd = paste(gw_filter_sh, rsid_file, gwas_file),
                  col.names = c("rsid", "beta", "se", "pval", "eaf", "ea", "oa"),
                  sep = "\t")
gwas_raw <- gwas_raw[!is.na(beta) & !is.na(se) & se > 0 & !is.na(rsid) & rsid != ""]
setkey(gwas_raw, rsid)
gwas_out <- data.frame(SNP = gwas_raw$rsid, beta.outcome = gwas_raw$beta,
                       se.outcome = gwas_raw$se, pval.outcome = gwas_raw$pval,
                       eaf.outcome = gwas_raw$eaf,
                       effect_allele.outcome = gwas_raw$ea,
                       other_allele.outcome = gwas_raw$oa,
                       outcome = "GCST90468131", id.outcome = "GCST90468131",
                       stringsAsFactors = FALSE)
# Phase 3 -- parallel PLINK clump + MR per metabolite.
met_ids <- unique(all_gws$met_id)
# Representative per-metabolite routine: clump GWS instruments, harmonize with
# the AR GWAS subset, and run IVW / weighted median / weighted mode / MR-RAPS
# (Wald ratio when only one instrument survives).
run_mr <- function(mid) {
  gws <- all_gws[met_id == mid]
  if (nrow(gws) == 0) return(NULL)
  exp_df <- data.frame(SNP = gws$rsid, beta.exposure = gws$beta,
                       se.exposure = gws$standard_error, pval.exposure = gws$p_value,
                       eaf.exposure = gws$effect_allele_frequency,
                       effect_allele.exposure = gws$effect_allele,
                       other_allele.exposure = gws$other_allele,
                       exposure = mid, id.exposure = mid, stringsAsFactors = FALSE)
  exp_clumped <- tryCatch(clump_data(exp_df, clump_r2 = 0.001, clump_kb = 500,
                                     pop = "EUR", bfile = bfile, plink_bin = plink_bin),
                          error = function(e) NULL)
  if (is.null(exp_clumped) || nrow(exp_clumped) == 0) return(NULL)
  gwas_sub <- gwas_out[gwas_out$SNP %in% exp_clumped$SNP, ]
  if (nrow(gwas_sub) == 0) return(NULL)
  dat <- harmonise_data(exp_clumped, gwas_sub, action = 2)
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) == 0) return(NULL)
  if (nrow(dat) == 1) {
    wr <- mr(dat, method_list = "mr_wald_ratio")
    return(data.frame(metabolite = mid, nsnp = 1, method = "Wald ratio",
                      beta = wr$b, se = wr$se, pval = wr$pval))
  }
  mr_res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  raps <- tryCatch(mr.raps(dat$beta.exposure, dat$beta.outcome,
                           dat$se.exposure, dat$se.outcome,
                           over.dispersion = TRUE, loss.function = "huber"),
                   error = function(e) NULL)
  rows <- data.frame(metabolite = mid, nsnp = mr_res$nsnp, method = mr_res$method,
                     beta = mr_res$b, se = mr_res$se, pval = mr_res$pval)
  if (!is.null(raps))
    rows <- rbind(rows, data.frame(metabolite = mid, nsnp = nrow(dat),
                                   method = "MR-RAPS", beta = raps$beta.hat,
                                   se = raps$beta.se,
                                   pval = 2 * pnorm(-abs(raps$beta.hat / raps$beta.se))))
  rows
}
results <- rbindlist(Filter(Negate(is.null),
                            mclapply(met_ids, run_mr, mc.cores = n_cores)), fill = TRUE)
# (elided: repeated clump+MR invocations across the significant metabolites; one representative routine is shown above)
ivw <- results[results$method %in% c("Inverse variance weighted", "Wald ratio"), ]
ivw$fdr <- p.adjust(ivw$pval, method = "BH")
# Expected: only GCST90199762 (N6-carbamoylthreonyladenosine) remains
# significant after multiple-testing correction (FDR < 0.05)
fwrite(results, file.path(res_dir, "C02_metabolite_AR_MR_all.tsv"), sep = "\t")
fwrite(ivw,     file.path(res_dir, "C02_metabolite_AR_MR_ivw.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 3. Mediation analysis (paper Figure 5c)
#    Indirect effect = beta_gm x beta_ma with delta-method SE; the total
#    gene -> AR estimate is the denominator. Expected: 6 significant chains
#    (FDR < 0.05), all through GCST90199762.
# -----------------------------------------------------------------------------
gene_met <- fread(file.path(res_dir, "C01_gene_metabolite_MR_ivw.tsv"))
met_ar   <- fread(file.path(res_dir, "C02_metabolite_AR_MR_ivw.tsv"))
gene_ar  <- fread(file.path(res_dir, "05_MR_ivw_fdr.tsv"))
gene_met_sig <- gene_met[fdr < 0.05]
met_ar_sig   <- met_ar[fdr < 0.05]
gene_ar_sig  <- gene_ar[fdr < 0.05 & method %in% c("Inverse variance weighted", "Wald ratio")]
# Best gene -> AR estimate per gene (lowest-p-value profile)
gene_ar_best <- gene_ar_sig[, .SD[which.min(pval)], by = gene][
  , .(gene, beta_total = beta, se_total = se)]
chains <- merge(gene_met_sig[, .(gene, metabolite, beta_gm = beta, se_gm = se)],
                met_ar_sig[,  .(metabolite, beta_ma = beta, se_ma = se)],
                by = "metabolite")
chains <- merge(chains, gene_ar_best, by = "gene")
if (nrow(chains) == 0) stop("No complete mediation chains")
chains[, beta_indirect := beta_gm * beta_ma]
chains[, se_indirect   := abs(beta_indirect) * sqrt((se_gm / beta_gm)^2 +
                                                     (se_ma / beta_ma)^2)]
chains[, z_indirect    := beta_indirect / se_indirect]
chains[, pval_indirect := 2 * pnorm(-abs(z_indirect))]
chains[, proportion_mediated := beta_indirect / beta_total]
chains[, fdr_indirect  := p.adjust(pval_indirect, method = "BH")]
sig_chains <- chains[fdr_indirect < 0.05]
fwrite(chains,     file.path(res_dir, "C03_mediation_all_chains.tsv"),  sep = "\t")
fwrite(sig_chains, file.path(res_dir, "C03_mediation_sig_chains.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 4. Colocalization of the mediating metabolite with the AR locus (paper
#    Figure 5d). The metabolite GWAS is a quantitative trait, so MAF and N
#    must be supplied in addition to beta/varbeta; the AR GWAS is case-control
#    (N = 394,626, s = 0.20). Shared-variant support: PP.H4/(PP.H3+PP.H4) > 0.7.
# -----------------------------------------------------------------------------
sig_chains <- fread(file.path(res_dir, "C03_mediation_sig_chains.tsv"))
if (nrow(sig_chains) == 0) stop("No significant chains to colocalize")
sig_mets <- unique(sig_chains$metabolite)  # GCST90199762
# AR GWAS chr/pos/rsid map for position-based joining with the metabolite
# (awk pre-filter at shell level; duplicate rsids from multi-allelic sites are
# reduced to the most significant |z| variant)
gwas_raw <- fread(cmd = paste0(
    "gunzip -c ", gwas_file,
    " | awk -F'\\t' 'NR==1{print \"chr\\tpos\\trsid\\tbeta\\tstandard_error\\teaf\"; next}",
    " $10!=\"\"&&$5!=\"\"&&$6!=\"\"{print $1\"\\t\"$2\"\\t\"$10\"\\t\"$5\"\\t\"$6\"\\t\"$7}'"),
  sep = "\t")
setnames(gwas_raw, "eaf", "effect_allele_frequency")
gwas_raw[, chr := as.character(chr)]
gwas_raw[, pos := as.integer(pos)]
gwas_raw <- gwas_raw[!is.na(beta) & !is.na(standard_error) & standard_error > 0]
gwas_raw[, varbeta := standard_error^2]
if (anyDuplicated(gwas_raw$rsid)) {
  gwas_raw[, absz := abs(beta) / standard_error]
  setorder(gwas_raw, -absz)
  gwas_raw <- gwas_raw[!duplicated(rsid)]
  gwas_raw[, absz := NULL]
}
setkey(gwas_raw, rsid)
N_gwas <- 394626L
s_gwas <- 0.20
coloc_results <- list()
for (mid in sig_mets) {
  mf <- file.path(met_dir, paste0(mid, "_buildGRCh38.tsv.gz"))
  if (!file.exists(mf)) next
  met <- tryCatch(fread(mf), error = function(e) NULL)
  if (is.null(met) || nrow(met) == 0) next
  # Map positions to rsids via the AR GWAS map (avoids a separate snp150 scan)
  if (!"rsid" %in% names(met)) {
    met[, chr := as.character(chromosome)]
    met[, pos := as.integer(base_pair_location)]
    met <- merge(met, gwas_raw[, .(chr, pos, rsid)],
                 by = c("chr", "pos"), all.x = FALSE)
  }
  met <- met[!is.na(rsid) & !is.na(beta) & !is.na(standard_error) & standard_error > 0]
  # (elided: duplicate-rsid deduplication; identical to the AR GWAS dedup above)
  if (nrow(met) == 0) next
  gwas_sub    <- gwas_raw[J(met$rsid), nomatch = 0]
  common_snps <- intersect(met$rsid, gwas_sub$rsid)
  if (length(common_snps) < 5) next
  met_sub  <- met[rsid %in% common_snps]
  gwas_sub <- gwas_sub[rsid %in% common_snps]
  coloc_res <- tryCatch(suppressMessages(coloc.abf(
    list(beta = met_sub$beta, varbeta = met_sub$standard_error^2,
         snp = met_sub$rsid, type = "quant", N = 8000L,
         MAF = met_sub$effect_allele_frequency),
    list(beta = gwas_sub$beta, varbeta = gwas_sub$varbeta,
         snp = gwas_sub$rsid, type = "cc", N = N_gwas, s = s_gwas)))$summary,
    error = function(e) NULL)
  if (is.null(coloc_res)) next
  coloc_results[[mid]] <- data.table(metabolite = mid,
    PP.H0 = unname(coloc_res["PP.H0.abf"]), PP.H1 = unname(coloc_res["PP.H1.abf"]),
    PP.H2 = unname(coloc_res["PP.H2.abf"]), PP.H3 = unname(coloc_res["PP.H3.abf"]),
    PP.H4 = unname(coloc_res["PP.H4.abf"]))
}
coloc_out <- rbindlist(coloc_results)
if (nrow(coloc_out) == 0) stop("No colocalization results")
coloc_out[, coloc_sig := PP.H4 / (PP.H3 + PP.H4) > 0.7]
# Adjudication: the posterior is dominated by H3 (PP.H3 = 0.984, PP.H4 = 0.003;
# ratio 0.003) -- distinct causal variants, so the shared-variant test is failed.
final_chains <- merge(sig_chains, coloc_out[, .(metabolite, PP.H4, coloc_sig)],
                      by = "metabolite")
fwrite(coloc_out,    file.path(res_dir, "C04_metabolite_AR_coloc.tsv"),    sep = "\t")
fwrite(final_chains, file.path(res_dir, "C04_final_mediation_chains.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 5. Paper Figure 5 panels
#    a gene-to-metabolite summary, b metabolite->AR volcano, c mediation
#    chains, d metabolite-AR coloc posteriors. PDF output via ggsave with the
#    high-contrast NPG palette.
# -----------------------------------------------------------------------------
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
col_risk    <- "#D62728"   # NPG red (risk-increasing)
col_protect <- "#1F77B4"   # NPG blue (protective)
col_core    <- "#FF7F0E"   # NPG orange (highlight)
col_sig     <- "#2CA02C"   # NPG green
col_grey    <- "#7F7F7F"
# Gene symbols for labeling
sym_dt <- fread(file.path(res_dir, "gene_symbol_full.tsv"))
setnames(sym_dt, c("gene", "gene_symbol"))
sym_lup <- setNames(sym_dt$gene_symbol, sym_dt$gene)
# --- Panel 5a: gene-to-metabolite MR summary ------------------------------
mediation <- fread(file.path(res_dir, "C03_mediation_sig_chains.tsv"))
c01 <- fread(file.path(res_dir, "C01_gene_metabolite_MR_ivw.tsv"))
c01[sym_dt, on = "gene", symbol := i.gene_symbol]
c01[, is_mediation := gene %in% unique(mediation$gene)]
c01[, is_sig := fdr < 0.05]
c01[, neg_log10_pval := -log10(pval)]
c01[, label := ifelse(is_mediation, symbol, NA_character_)]
p_a <- ggplot(c01, aes(x = beta, y = neg_log10_pval)) +
  geom_point(data = c01[is_sig == FALSE], aes(color = "Non-significant"),
             size = 0.3, alpha = 0.4) +
  geom_point(data = c01[is_sig == TRUE & is_mediation == FALSE],
             aes(color = "Significant"), size = 0.5, alpha = 0.7) +
  geom_point(data = c01[is_mediation == TRUE],
             aes(color = "Mediation chain"), size = 2, alpha = 1) +
  geom_text_repel(data = c01[is_mediation == TRUE], aes(label = label),
                  size = 2.5, fontface = "italic", box.padding = 0.5,
                  max.overlaps = 20, min.segment.length = 0.2, color = col_core) +
  scale_color_manual(name = NULL,
    values = c("Non-significant" = col_grey, "Significant" = col_sig,
               "Mediation chain" = col_core),
    breaks = c("Mediation chain", "Significant", "Non-significant")) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "grey50") +
  labs(title = "Gene -> Metabolite MR results",
       subtitle = "1,693 significant pairs | 126 genes -> 314 metabolites (FDR < 0.05)",
       x = "MR effect on metabolite level (beta, IVW)",
       y = expression(-log[10](italic(P)))) +
  theme_nature() +
  theme(plot.subtitle = element_text(size = 7, color = "grey40"),
        legend.position = "bottom", legend.direction = "horizontal")
ggsave(file.path(fig_dir, "Fig5a_gene_metabolite_MR.pdf"), p_a,
       width = 6, height = 5, device = "pdf")
# --- Panel 5b: metabolite-to-AR MR volcano ---------------------------------
metab_mr <- fread(file.path(res_dir, "C02_metabolite_AR_MR_ivw.tsv"))
metab_mr[, log10fdr := -log10(fdr)]
p_b <- ggplot(metab_mr[order(fdr)], aes(x = beta, y = log10fdr)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = col_risk, linewidth = 0.5) +
  geom_point(data = metab_mr[fdr >= 0.05], aes(size = nsnp),
             color = "grey60", alpha = 0.5, shape = 16) +
  geom_point(data = metab_mr[fdr < 0.05], aes(size = nsnp),
             color = col_risk, alpha = 0.95, shape = 16) +
  geom_text_repel(data = metab_mr[fdr < 0.05], aes(label = metabolite),
                  size = 2.8, max.overlaps = 10, fontface = "bold",
                  segment.size = 0.3, segment.color = "grey50",
                  box.padding = 0.5, point.padding = 0.3) +
  scale_size_continuous(range = c(2, 6), name = "Number of\ninstruments",
                        breaks = c(2, 5, 10, 20, 50, 100)) +
  labs(title = "Plasma metabolite MR effects on allergic rhinitis",
       subtitle = sprintf("%d metabolites tested | 1 significant (FDR < 0.05): GCST90199762",
                          nrow(metab_mr)),
       x = expression("MR effect on allergic rhinitis (" * beta * ")"),
       y = expression(-log[10]("FDR"))) +
  theme_nature()
ggsave(file.path(fig_dir, "Fig5b_metabolite_AR_volcano.pdf"), p_b,
       width = 6.5, height = 5, device = "pdf")
# --- Panel 5c: six mediation chains through GCST90199762 --------------------
med_sig <- copy(mediation)
med_sig[, gene_symbol := coalesce(sym_lup[gene], gene)]
med_sig[, prop_pct := proportion_mediated * 100]
med_sig[, chain_label := paste0(gene_symbol, "  ->  GCST90199762")]
med_sig <- med_sig[order(beta_indirect)]
med_sig[, chain_f := fct_reorder(chain_label, beta_indirect)]
p_c1 <- ggplot(med_sig, aes(x = prop_pct, y = chain_f)) +
  geom_segment(aes(x = 0, xend = prop_pct, yend = chain_f),
               color = "grey60", linewidth = 0.8) +
  geom_point(aes(color = prop_pct), size = 3.5) +
  scale_color_gradient(low = "#AEC6E8", high = col_risk, name = "% mediated") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Proportion of total effect mediated via metabolite",
       x = "Proportion mediated (%)", y = NULL) +
  theme_nature() + theme(axis.text.y = element_text(face = "italic", size = 8))
p_c2 <- ggplot(med_sig, aes(y = gene_symbol)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = beta_total, yend = gene_symbol),
               color = "grey70", linewidth = 1.5, alpha = 0.4) +
  geom_point(aes(x = beta_total, color = "Total"), size = 3) +
  geom_point(aes(x = beta_indirect, color = "Indirect"), size = 2.5, shape = 17) +
  scale_color_manual(values = c("Total" = col_protect, "Indirect" = col_risk),
                     name = "Effect") +
  labs(title = "Total vs. indirect effect estimates",
       x = expression("MR effect size (" * beta * ")"), y = NULL) +
  theme_nature() + theme(axis.text.y = element_text(face = "italic", size = 8))
p_c <- p_c1 | p_c2
p_c <- p_c + plot_annotation(
  title = "Gene -> Metabolite -> Allergic rhinitis mediation chains",
  subtitle = sprintf("n = %d significant indirect pathways (FDR < 0.05) | all via GCST90199762",
                     nrow(med_sig)),
  theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)))
ggsave(file.path(fig_dir, "Fig5c_mediation_chains.pdf"), p_c,
       width = 12, height = 5, device = "pdf")
# --- Panel 5d: metabolite-AR colocalization posteriors ----------------------
# Posterior probabilities of GCST90199762 x AR GWAS colocalization; H3
# (independent signals) dominates over H4 (shared causal variant).
pp_vals <- c(0.000, 0.000, 0.013, 0.984, 0.003)
pp_df <- data.frame(Hypothesis = factor(paste0("H", 0:4), levels = paste0("H", 0:4)),
                    PP = pp_vals,
                    label = c("H0: No association", "H1: Metabolite only",
                              "H2: AR only", "H3: Independent signals",
                              "H4: Shared causal variant"),
                    stringsAsFactors = FALSE)
pp_df$highlight <- pp_df$Hypothesis == "H3"
p_d1 <- ggplot(pp_df, aes(x = Hypothesis, y = PP, fill = highlight)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", PP)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("FALSE" = "#AEC6E8", "TRUE" = col_core), guide = "none") +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  labs(title = "Metabolite-AR colocalization: posterior probabilities",
       subtitle = "GCST90199762 (plasma metabolite) x Allergic rhinitis GWAS",
       x = "Bayesian hypothesis", y = "Posterior probability") +
  theme_nature()
# Pathway schematic with the two-step MR mediation summary
p_d2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.8, label = "Gene expression\n(CD4+ T cells)",
           size = 3.5, fontface = "bold", hjust = 0.5) +
  annotate("segment", x = 0.6, xend = 1.2, y = 0.8, yend = 0.8,
           arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.8) +
  annotate("text", x = 1.35, y = 0.8, label = "Plasma\nMetabolite\n(GCST90199762)",
           size = 3.5, fontface = "bold", hjust = 0.5) +
  annotate("segment", x = 1.5, xend = 2.1, y = 0.8, yend = 0.8,
           arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.8,
           linetype = "dashed", color = "grey50") +
  annotate("text", x = 2.25, y = 0.8, label = "Allergic\nRhinitis",
           size = 3.5, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 0.9, y = 0.6,
           label = paste0("6 gene chains\nFDR < 0.05\nmediation proportion: 8-33%"),
           size = 2.5, color = "grey40") +
  annotate("text", x = 1.8, y = 0.6,
           label = expression(beta * " = -0.211\nFDR = 5.84e-04\n(metabolite -> AR)"),
           size = 2.5, color = "grey40") +
  annotate("text", x = 1.35, y = 0.35,
           label = "Coloc PP.H3 = 0.984\n(independent signals -- mediation\nnot supported by coloc)",
           size = 2.5, color = col_risk, hjust = 0.5, fontface = "italic") +
  xlim(0, 2.75) + ylim(0.1, 1.05) +
  labs(title = "Two-step MR mediation pathway") +
  theme_void() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
p_d <- p_d1 | p_d2
ggsave(file.path(fig_dir, "Fig5d_metabolite_coloc.pdf"), p_d,
       width = 11, height = 4.5, device = "pdf")
