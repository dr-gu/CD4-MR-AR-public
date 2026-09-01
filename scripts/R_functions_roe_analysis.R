# =============================================================================
# R_functions_roe_analysis.R
# Shared helper library for the single-cell module (Results 2.5):
# Ro/e positive-enrichment analysis, Fisher tests, miloR neighborhood
# differential abundance, and figure themes (theme_nature, NPG palette).
# Paper: Genetic prioritization of CD4+ T cell-expressed candidate genes for
#        allergic rhinitis (Biomedicines)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(stringr)
})

# Global constants

NPG_RED    <- "#D62728"
NPG_BLUE   <- "#1F77B4"
NPG_GREEN  <- "#2CA02C"
NPG_ORANGE <- "#FF7F0E"
NPG_PURPLE <- "#9467BD"

# Nature theme (consistent with sc_10/sc_15/sc_21)

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

# =============================================================================
# 1. Data loading
# =============================================================================

#' Load target gene list (with tier classification and MR direction)
#'
#' @param base_dir Project root directory (repo root). The default assumes the
#'   working directory is the repo root.
#' @return data.table: gene_symbol, tier, direction, beta, fdr
load_target_genes <- function(base_dir = ".") {
  # 129 significant genes
  sig_genes <- fread(file.path(base_dir, "results", "mr_sig_gene_symbols.txt"),
                     header = FALSE, col.names = "gene_symbol")

  # Core genes + MR parameters
  core_dt <- fread(file.path(base_dir, "tables", "Table1_core10_genes_AR.tsv"))

  # Full MR results (beta, fdr)
  mr_full <- fread(file.path(base_dir, "results", "05_MR_ivw_fdr.tsv"))
  sym_dt  <- fread(file.path(base_dir, "results", "gene_symbol_full.tsv"))
  setnames(sym_dt, c("gene", "gene_symbol"))
  sym_lup <- setNames(sym_dt$gene_symbol, sym_dt$gene)

  mr_best <- mr_full[fdr < 0.05][order(fdr, pval)][!duplicated(gene)]
  mr_best[, gene_symbol := coalesce(sym_lup[gene], gene)]
  mr_best[, direction := ifelse(beta > 0, "Risk", "Protect")]

  # Tier stratification (based on MR_04_downstream_research_v3)
  tier1 <- c("TSLP", "ORMDL3", "SLC25A46", "TLR1")
  tier2 <- c("IL18R1", "IL18RAP", "IKZF3", "CD247")
  tier3 <- c("PDCD1", "D2HGDH", "STAT6", "LMBRD2", "IL1R2")

  # Merge
  gene_dt <- data.table(gene_symbol = sig_genes$gene_symbol)
  gene_dt[, tier := fifelse(gene_symbol %in% tier1, "Tier1",
                     fifelse(gene_symbol %in% tier2, "Tier2",
                     fifelse(gene_symbol %in% tier3, "Tier3", "Other")))]
  gene_dt[, core := gene_symbol %in% core_dt$Gene]

  # Merge direction
  gene_dt <- merge(gene_dt,
                   mr_best[, .(gene_symbol, direction, beta, fdr)],
                   by = "gene_symbol", all.x = TRUE)
  gene_dt[is.na(direction), direction := "Unknown"]

  # Sort: tier -> direction -> gene name
  gene_dt[, tier_order := fifelse(tier == "Tier1", 1,
                           fifelse(tier == "Tier2", 2,
                           fifelse(tier == "Tier3", 3, 4)))]
  setorder(gene_dt, tier_order, gene_symbol)

  return(gene_dt[])
}

# =============================================================================
# 2. Ro/e calculation (adapted from tutu distribution_Roe.R)
# =============================================================================

#' Single-gene Ro/e calculation - cell type x condition distribution
#' preference of positive cells
#'
#' Adapted from the chi-squared observed/expected method in tutu
#' distribution_Roe.R. Two input modes are supported:
#'   1) a pre-extracted data.frame (expr_df) with the gene column plus
#'      celltype_col and condition_col
#'   2) a Seurat object (sce), extracted via FetchData (for quick testing only)
#'
#' @param sce optional Seurat object
#' @param expr_df optional pre-extracted data.frame (preferred, avoids passing
#'   large objects to parallel workers)
#' @param gene gene symbol
#' @param celltype_col cell type column name
#' @param condition_col condition column name
#' @param assay Seurat assay (sce mode only)
#' @param layer Seurat layer (sce mode only)
#' @return list(roe_matrix, observed, expected, n_pos, pct_pos)
calc_gene_roe <- function(sce = NULL,
                          expr_df = NULL,
                          gene,
                          celltype_col = "cell_type",
                          condition_col = "condition",
                          assay = "RNA",
                          layer = "data") {

  # Prefer the pre-extracted data.frame
  if (!is.null(expr_df)) {
    if (!gene %in% names(expr_df)) {
      return(list(roe = NULL, observed = NULL, expected = NULL,
                  n_pos = 0, pct_pos = 0, error = "Gene not in expr_df"))
    }
    expr_data <- expr_df[, c(gene, celltype_col, condition_col)]
    names(expr_data)[1] <- "expression"
  } else if (!is.null(sce)) {
    if (!gene %in% rownames(sce)) {
      return(list(roe = NULL, observed = NULL, expected = NULL,
                  n_pos = 0, pct_pos = 0, error = "Gene not in Seurat object"))
    }
    expr_data <- tryCatch(
      FetchData(sce, vars = c(gene, celltype_col, condition_col),
                assay = assay, layer = layer),
      error = function(e) NULL
    )
    if (is.null(expr_data)) {
      return(list(roe = NULL, observed = NULL, expected = NULL,
                  n_pos = 0, pct_pos = 0, error = "FetchData failed"))
    }
    names(expr_data)[1] <- "expression"
  } else {
    return(list(roe = NULL, observed = NULL, expected = NULL,
                n_pos = 0, pct_pos = 0, error = "No data provided"))
  }

  # Binarize
  expr_data$positive <- expr_data$expression > 0
  n_total <- nrow(expr_data)
  n_pos <- sum(expr_data$positive)
  pct_pos <- mean(expr_data$positive) * 100

  if (n_pos < 5) {
    return(list(roe = NULL, observed = NULL, expected = NULL,
                n_pos = n_pos, pct_pos = pct_pos,
                error = "Too few positive cells (< 5)"))
  }

  # Keep only positive cells
  pos_cells <- expr_data[expr_data$positive, ]

  # Build contingency table
  pos_cells[[celltype_col]] <- as.factor(pos_cells[[celltype_col]])
  pos_cells[[condition_col]] <- as.factor(pos_cells[[condition_col]])

  cont_tab <- tryCatch(
    xtabs(as.formula(paste0("~", celltype_col, "+", condition_col)),
          data = pos_cells),
    error = function(e) NULL
  )
  if (is.null(cont_tab) || nrow(cont_tab) < 2) {
    return(list(roe = NULL, observed = NULL, expected = NULL,
                n_pos = n_pos, pct_pos = pct_pos,
                error = "Contingency table too small"))
  }

  # Chi-squared test to extract observed/expected
  chisq <- tryCatch(
    chisq.test(as.matrix(cont_tab)),
    error = function(e) NULL
  )
  if (is.null(chisq)) {
    observed <- as.matrix(cont_tab)
    rs <- rowSums(observed)
    cs <- colSums(observed)
    expected <- outer(rs, cs) / sum(observed)
    dimnames(expected) <- dimnames(observed)
    roe <- observed / expected
  } else {
    roe <- chisq$observed / chisq$expected
    observed <- chisq$observed
    expected <- chisq$expected
  }

  return(list(
    roe        = roe,
    observed   = observed,
    expected   = expected,
    n_pos      = n_pos,
    pct_pos    = pct_pos,
    n_total    = n_total,
    error      = NULL
  ))
}

#' Convert a Ro/e matrix to long-format data.table
#'
#' @param roe_result return value of calc_gene_roe
#' @param gene gene symbol
#' @return data.table: gene, cell_type, condition, observed, expected, roe
roe_matrix_to_dt <- function(roe_result, gene) {
  if (!is.null(roe_result$error) || is.null(roe_result$roe)) {
    return(data.table(
      gene = gene, cell_type = NA_character_, condition = NA_character_,
      observed = NA_real_, expected = NA_real_, roe = NA_real_,
      error_msg = roe_result$error %||% "unknown"
    ))
  }

  roe_mat <- roe_result$roe
  obs_mat <- roe_result$observed
  exp_mat <- roe_result$expected

  cell_types <- rownames(roe_mat)
  conditions  <- colnames(roe_mat)

  result <- CJ(cell_type = cell_types, condition = conditions)
  result[, gene := gene]

  for (i in seq_len(nrow(result))) {
    ct <- result$cell_type[i]
    cd <- result$condition[i]
    result[i, `:=`(
      observed = obs_mat[ct, cd],
      expected = exp_mat[ct, cd],
      roe      = roe_mat[ct, cd]
    )]
  }

  if ("error_msg" %in% names(result)) result[, error_msg := NULL]
  return(result[])
}

# =============================================================================
# 3. Fisher exact test (layer 2 - comparing positive proportions within
#    each cell type)
# =============================================================================

#' Single-gene Fisher exact test across all cell types
#'
#' For each cell type, build a 2x2 contingency table
#' (positive/negative x condition1/condition2) and run a Fisher exact test
#' comparing positive cell proportions.
#'
#' @param sce Seurat object
#' @param gene gene symbol
#' @param celltype_col cell type column
#' @param condition_col condition column
#' @param min_pos_total minimum total positive cells (filter threshold)
#' @param min_pos_per_group minimum positive cells per group
#' @return data.table
calc_gene_fisher <- function(sce = NULL,
                             expr_df = NULL,
                             gene,
                             celltype_col = "cell_type",
                             condition_col = "condition",
                             min_pos_total = 5,
                             min_pos_per_group = 3) {

  # Prefer the pre-extracted data.frame
  if (!is.null(expr_df)) {
    if (!gene %in% names(expr_df)) {
      return(data.table(gene = gene, n_tested = 0, error = "Gene not in expr_df"))
    }
    expr_data <- expr_df[, c(gene, celltype_col, condition_col)]
    names(expr_data)[1] <- "expression"
  } else if (!is.null(sce)) {
    if (!gene %in% rownames(sce)) {
      return(data.table(gene = gene, n_tested = 0, error = "Gene not found"))
    }
    expr_data <- tryCatch(
      FetchData(sce, vars = c(gene, celltype_col, condition_col)),
      error = function(e) NULL
    )
    if (is.null(expr_data)) {
      return(data.table(gene = gene, n_tested = 0, error = "FetchData failed"))
    }
    names(expr_data)[1] <- "expression"
  } else {
    return(data.table(gene = gene, n_tested = 0, error = "No data provided"))
  }

  expr_data$positive <- expr_data$expression > 0
  conditions <- unique(expr_data[[condition_col]])
  if (length(conditions) != 2) {
    return(data.table(gene = gene, n_tested = 0,
                      error = "Need exactly 2 condition levels"))
  }
  cond1 <- conditions[1]
  cond2 <- conditions[2]

  cell_types <- unique(expr_data[[celltype_col]])
  result_list <- list()

  for (ct in cell_types) {
    ct_data <- expr_data[expr_data[[celltype_col]] == ct, ]

    d1 <- ct_data[ct_data[[condition_col]] == cond1, ]
    d2 <- ct_data[ct_data[[condition_col]] == cond2, ]

    n1_pos <- sum(d1$positive)
    n1_neg <- sum(!d1$positive)
    n2_pos <- sum(d2$positive)
    n2_neg <- sum(!d2$positive)

    n1_total <- nrow(d1)
    n2_total <- nrow(d2)
    total_pos <- n1_pos + n2_pos

    # Filter
    if (total_pos < min_pos_total ||
        n1_pos < min_pos_per_group ||
        n2_pos < min_pos_per_group) {
      next
    }

    # Proportions
    pct1 <- n1_pos / n1_total * 100
    pct2 <- n2_pos / n2_total * 100

    # Ro/e (within cell type)
    expected1 <- (n1_total / (n1_total + n2_total)) * total_pos
    expected2 <- (n2_total / (n1_total + n2_total)) * total_pos
    roe1 <- if (expected1 > 0) n1_pos / expected1 else NA
    roe2 <- if (expected2 > 0) n2_pos / expected2 else NA
    log2_roe_ratio <- log2((n1_pos + 1) / (expected1 + 1) /
                           ((n2_pos + 1) / (expected2 + 1)))

    # Fisher's exact test
    cont_table <- matrix(c(n1_pos, n2_pos, n1_neg, n2_neg), nrow = 2)
    ft <- tryCatch(
      fisher.test(cont_table, alternative = "two.sided"),
      error = function(e) NULL
    )

    result_list[[length(result_list) + 1]] <- data.table(
      gene          = gene,
      cell_type     = ct,
      n_cells       = n1_total + n2_total,
      n_cond1       = n1_total,
      n_cond2       = n2_total,
      n_cond1_pos   = n1_pos,
      n_cond1_neg   = n1_neg,
      n_cond2_pos   = n2_pos,
      n_cond2_neg   = n2_neg,
      pct_cond1_pos = round(pct1, 3),
      pct_cond2_pos = round(pct2, 3),
      roe_cond1     = round(roe1, 4),
      roe_cond2     = round(roe2, 4),
      log2_roe_ratio = round(log2_roe_ratio, 4),
      odds_ratio    = if (!is.null(ft)) round(ft$estimate, 4) else NA_real_,
      or_ci_lower   = if (!is.null(ft)) round(ft$conf.int[1], 4) else NA_real_,
      or_ci_upper   = if (!is.null(ft)) round(ft$conf.int[2], 4) else NA_real_,
      fisher_pval   = if (!is.null(ft)) ft$p.value else NA_real_
    )
  }

  if (length(result_list) == 0) {
    return(data.table(gene = gene, n_tested = 0))
  }
  return(rbindlist(result_list))
}

# =============================================================================
# 4. Visualization - Ro/e heatmap (adapted from the ggplot pattern in tutu
#    distribution_Roe)
# =============================================================================

#' Multi-gene Ro/e heatmap panel (tutu style)
#'
#' For each gene, draw a cell_type x condition Ro/e heatmap using an NPG
#' red-white-blue gradient. Multi-gene panels are arranged with patchwork.
#'
#' @param roe_list named list of calc_gene_roe results
#' @param max_roe Ro/e cap (clipping)
#' @param ncol number of panel columns
#' @param tile_color tile border color
#' @return ggplot object
plot_roe_heatmap <- function(roe_list,
                             max_roe = 3,
                             ncol = 4,
                             tile_color = "grey90") {

  plot_list <- list()
  gene_names <- names(roe_list)

  for (gn in gene_names) {
    res <- roe_list[[gn]]
    if (!is.null(res$error) || is.null(res$roe)) next

    dt <- roe_matrix_to_dt(res, gn)
    if (nrow(dt) == 0 || all(is.na(dt$roe))) next

    # Clip
    dt[, roe_clipped := pmax(pmin(roe, max_roe), 1 / max_roe)]
    dt[, roe_clipped := fifelse(is.na(roe_clipped), 1, roe_clipped)]

    n_pos <- res$n_pos
    pct   <- round(res$pct_pos, 1)

    p <- ggplot(dt, aes(x = condition, y = cell_type)) +
      geom_tile(aes(fill = roe_clipped), color = tile_color) +
      scale_fill_gradient2(
        name = "Ro/e",
        low = NPG_BLUE, mid = "white", high = NPG_RED,
        midpoint = 1,
        limits = c(1 / max_roe, max_roe),
        oob = scales::squish
      ) +
      scale_y_discrete(expand = c(0, 0), position = "right") +
      scale_x_discrete(expand = c(0, 0)) +
      labs(
        title = gn,
        subtitle = sprintf("Pos: %d (%.1f%%)", n_pos, pct)
      ) +
      theme_nature(9) +
      theme(
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y.right = element_text(size = 7, color = "black"),
        axis.text.x.bottom = element_text(size = 7, color = "black",
                                          angle = 45, hjust = 1),
        legend.position = "right",
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7, color = "grey50")
      )

    plot_list[[gn]] <- p
  }

  if (length(plot_list) == 0) return(NULL)
  wrap_plots(plot_list, ncol = ncol) +
    plot_annotation(theme = theme(plot.title = element_text(size = 12, face = "bold")))
}

#' Bubble plot of significant Fisher results
#'
#' @param fisher_dt merged Fisher result data.table (requires fdr column)
#' @param fdr_threshold significance threshold
#' @param max_genes maximum number of genes to show
#' @return ggplot object
plot_fisher_bubble <- function(fisher_dt, fdr_threshold = 0.05, max_genes = 40) {

  sig_dt <- fisher_dt[fdr < fdr_threshold]
  if (nrow(sig_dt) == 0) {
    message("No significant gene-cell_type pairs at FDR < ", fdr_threshold)
    return(NULL)
  }

  # Select top genes
  top_genes <- sig_dt[, .(n_sig = .N), by = gene][order(-n_sig)]
  if (nrow(top_genes) > max_genes) top_genes <- top_genes[1:max_genes]

  plot_dt <- sig_dt[gene %in% top_genes$gene]
  plot_dt[, neg_log10_fdr := -log10(fdr)]
  plot_dt[, enrich_group := fifelse(log2_roe_ratio > 0, "Pre", "Post")]

  # Sort by hierarchical clustering
  gene_order <- plot_dt[, .(mean_roe = mean(abs(log2_roe_ratio))), by = gene][order(-mean_roe)]$gene
  ct_order   <- plot_dt[, .(mean_roe = mean(abs(log2_roe_ratio))), by = cell_type][order(-mean_roe)]$cell_type

  plot_dt[, gene := factor(gene, levels = gene_order)]
  plot_dt[, cell_type := factor(cell_type, levels = ct_order)]

  ggplot(plot_dt, aes(x = gene, y = cell_type)) +
    geom_point(aes(size = neg_log10_fdr, color = log2_roe_ratio)) +
    scale_color_gradient2(
      name = "log2(Ro/e ratio)",
      low = NPG_BLUE, mid = "white", high = NPG_RED,
      midpoint = 0
    ) +
    scale_size_continuous(name = "-log10(FDR)", range = c(1, 6)) +
    labs(x = "", y = "", title = "Gene positive cell enrichment (FDR < 0.05)") +
    theme_nature(9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    )
}

#' Single-gene Ro/e detail plot (bar chart)
#'
#' @param roe_dt output of roe_matrix_to_dt
#' @param gene gene name
#' @return ggplot object
plot_gene_roe_detail <- function(roe_dt, gene) {

  dt <- copy(roe_dt)
  dt[, roe_label := round(roe, 2)]
  dt[is.na(roe_label), roe_label := 0]

  # log2(Ro/e) is more intuitive
  dt[, log2_roe := log2(pmax(roe, 0.125))]

  ggplot(dt, aes(x = cell_type, y = log2_roe, fill = condition)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("Pre" = NPG_RED, "Post" = NPG_BLUE)) +
    labs(
      title = gene,
      x = "", y = "log2(Ro/e)",
      subtitle = "Ro/e > 0 = enriched (observed > expected)"
    ) +
    theme_nature(9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top"
    )
}

# =============================================================================
# 5. miloR neighborhood differential abundance (layer 3 - based on horse
#    day2_3_milor.R)
# =============================================================================

#' Run miloR differential abundance analysis
#'
#' Standard miloR workflow (based on horse day2_3_milor.R):
#' Seurat->SCE->Milo->buildGraph->makeNhoods->countCells->calcNhoodDistance->
#' testNhoods->annotateNhoods
#'
#' @param sce Seurat object
#' @param condition_col condition column name
#' @param sample_col sample column name (used for countCells)
#' @param reduced_dim reduced dimension name (default "HARMONY")
#' @param k KNN parameter
#' @param d number of dimensions
#' @param prop neighborhood proportion
#' @return list(milo, da_results)
run_milor_da <- function(sce,
                          condition_col = "condition",
                          sample_col = "patient",
                          reduced_dim = "HARMONY",
                          k = 30, d = 20, prop = 0.2) {

  # Check packages
  if (!requireNamespace("miloR", quietly = TRUE)) {
    stop("miloR not installed. Run: BiocManager::install('miloR')")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment not installed.")
  }
  if (!requireNamespace("scater", quietly = TRUE)) {
    stop("scater not installed. Run: BiocManager::install('scater')")
  }

  library(miloR)
  library(SingleCellExperiment)
  library(scater)
  library(Matrix)

  # Matrix compatibility patch (see horse day2_3_milor.R)
  setAs("dgTMatrix", "dgCMatrix", function(from) {
    Matrix::sparseMatrix(
      i = from@i + 1L,
      j = from@j + 1L,
      x = from@x,
      dims = from@Dim,
      dimnames = from@Dimnames,
      giveCsparse = TRUE
    )
  })

  # Prepare metadata
  sce$orig.ident <- as.character(sce@meta.data[[sample_col]])
  sce$group <- factor(sce@meta.data[[condition_col]])

  # Seurat → SingleCellExperiment
  cat("  Converting to SingleCellExperiment...\n")
  sce_milo <- as.SingleCellExperiment(sce)

  # Use Harmony embedding as reducedDim
  if (reduced_dim %in% names(sce@reductions)) {
    reducedDim(sce_milo, "HARMONY") <- Embeddings(sce, reduced_dim)
  }

  # Milo workflow
  cat("  Building Milo object...\n")
  milo <- miloR::Milo(sce_milo)

  cat(sprintf("  Building KNN graph (k=%d, d=%d)...\n", k, d))
  milo <- miloR::buildGraph(milo, k = k, d = d)

  cat(sprintf("  Making neighborhoods (prop=%.1f, k=%d, d=%d)...\n", prop, k, d))
  milo <- makeNhoods(milo, prop = prop, k = k, d = d, refined = TRUE)

  cat("  Counting cells per neighborhood...\n")
  milo <- countCells(milo,
                     meta.data = data.frame(colData(milo)),
                     sample = "orig.ident")

  # Design matrix
  traj_design <- data.frame(colData(milo))[, c("orig.ident", "group")]
  traj_design$orig.ident <- as.factor(traj_design$orig.ident)
  traj_design <- distinct(traj_design)
  rownames(traj_design) <- traj_design$orig.ident

  cat("  Calculating neighborhood distances...\n")
  milo <- calcNhoodDistance(milo, d = d)

  cat("  Running differential abundance test...\n")
  da_results <- testNhoods(milo,
                           design = as.formula("~ group"),
                           design.df = traj_design)

  # Annotate neighborhoods
  cat("  Annotating neighborhoods...\n")
  da_results <- annotateNhoods(milo, da_results, coldata_col = "cell_type")

  cat(sprintf("  Done. %d neighborhoods, %d significant (SpatialFDR < 0.05)\n",
              nrow(da_results),
              sum(da_results$SpatialFDR < 0.05, na.rm = TRUE)))

  return(list(milo = milo, da_results = da_results))
}

#' miloR beeswarm plot (plotDAbeeswarm wrapper)
#'
#' @param da_results testNhoods output (after annotateNhoods)
#' @param group_by grouping column name
#' @param alpha transparency
#' @return ggplot object
plot_milor_beeswarm <- function(da_results,
                                 group_by = "cell_type",
                                 alpha = 0.8) {

  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    install.packages("ggbeeswarm", repos = "https://cran.rstudio.com")
  }
  library(ggbeeswarm)
  library(scales)

  plotDAbeeswarm(da_results, group.by = group_by, alpha = alpha) +
    scale_color_gradient2(
      midpoint = 0,
      low = NPG_BLUE, mid = "white", high = NPG_RED,
      space = "Lab"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "", y = "Log2 Fold Change") +
    theme_nature(9) +
    theme(axis.text = element_text(color = "black"))
}

#' miloR neighborhood graph plot
#'
#' @param milo Milo object
#' @param da_results DA results
#' @return ggplot object
plot_milor_nhood_graph <- function(milo, da_results) {

  if (!requireNamespace("miloR", quietly = TRUE)) {
    stop("miloR not installed")
  }
  library(miloR)
  library(scales)

  milo <- buildNhoodGraph(milo)
  plotNhoodGraphDA(milo, da_results) +
    scale_fill_gradient2(
      low = NPG_BLUE, mid = "lightgrey", high = NPG_RED,
      name = "log2FC",
      limits = c(-3, 3),
      oob = squish
    ) +
    theme_nature(9)
}

# =============================================================================
# 6. Helper functions
# =============================================================================

#' Print analysis summary
#'
#' @param roe_results named list of calc_gene_roe results
#' @param fisher_dt merged Fisher results
#' @param fdr_threshold significance threshold
print_analysis_summary <- function(roe_results, fisher_dt, fdr_threshold = 0.05) {
  n_genes <- length(roe_results)
  n_success <- sum(sapply(roe_results, function(x) is.null(x$error)))
  n_error  <- sum(sapply(roe_results, function(x) !is.null(x$error)))

  cat(sprintf("\n╔══════════════════════════════════════════╗\n"))
  cat(sprintf("║  Ro/e + Fisher Analysis Summary          ║\n"))
  cat(sprintf("╠══════════════════════════════════════════╣\n"))
  cat(sprintf("║  Genes analyzed:          %4d         ║\n", n_genes))
  cat(sprintf("║  Successful:              %4d         ║\n", n_success))
  cat(sprintf("║  Failed:                  %4d         ║\n", n_error))

  if (nrow(fisher_dt) > 0 && "fdr" %in% names(fisher_dt)) {
    n_sig <- sum(fisher_dt$fdr < fdr_threshold, na.rm = TRUE)
    n_genes_sig <- uniqueN(fisher_dt[fdr < fdr_threshold]$gene)
    n_ct_sig <- uniqueN(fisher_dt[fdr < fdr_threshold]$cell_type)
    cat(sprintf("║  Sig gene-CT pairs (FDR<%.2f): %4d   ║\n", fdr_threshold, n_sig))
    cat(sprintf("║  Sig unique genes:        %4d         ║\n", n_genes_sig))
    cat(sprintf("║  Sig unique cell types:   %4d         ║\n", n_ct_sig))
  }
  cat(sprintf("╚══════════════════════════════════════════╝\n\n"))
}

message("R_functions_roe_analysis.R loaded successfully.")
