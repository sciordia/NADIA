# =============================================================================
# Batch Correction — PVCA + BERT
# =============================================================================
#
# Two complementary modules for batch effect analysis:
#
# SECTION 1-2: PVCA (Principal Variance Component Analysis)
#   Decomposes total variance across experimental factors to identify the main
#   sources of variation (batch effects, biological covariates, etc.).
#   Public functions:
#     pvca_compute()               : Core variance decomposition (data.frame)
#     pvca_plot()                  : Bar plot of variance components (ggplot2)
#     pvca_analysis()              : Orchestrator (compute + plot + export)
#
# SECTION 3-4: Batch Correction (BERT)
#   Optional batch effect correction using BERT (ComBat/limma/ref with
#   hierarchical tree-based correction for missing-value-tolerant correction).
#   Public functions:
#     batch_correct_proteomics()   : Correct batch effects, add assay to SE
#
# SECTION 5: PCA Colored by Covariates
#   PCA scatter plots colored by different covariates for visual batch
#   effect identification.
#   Public functions:
#     pca_covariates_plot()        : Grid + individual PCA plots per covariate
#
# References:
#   PVCA: Li et al. (2009) Biostatistics 10(2):317-326
#         Cuklina et al. (2021) Molecular & Cellular Proteomics 20 (proBatch)
#   BERT: Habarta et al. (2025) Nature Communications 16:2044
#
# Implementation notes:
#   - PVCA uses lme4 directly (not pvca package) to allow automatic detection
#     and exclusion of problematic interaction terms where n_levels >= n_obs
#   - BERT handles missing values via hierarchical tree decomposition,
#     suitable for pre-imputation batch correction on protein-level data
#
# Dependencies:
#   PVCA: SummarizedExperiment, ggplot2, lme4
#   Batch Correction: SummarizedExperiment, BERT (Bioconductor)
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# =============================================================================
# SECTION 1: INTERNAL HELPERS
# =============================================================================

#' Validate that factors exist in colData
#' @keywords internal
.pvca_validate_factors <- function(se, factors) {
  cd_cols <- colnames(SummarizedExperiment::colData(se))
  missing <- setdiff(factors, cd_cols)
  if (length(missing) > 0)
    stop("Factor(s) not found in colData: ", paste(missing, collapse = ", "),
         "\n  Available columns: ", paste(cd_cols, collapse = ", "))
  invisible(TRUE)
}

#' Extract and prepare the abundance matrix, handling NAs
#' @keywords internal
.pvca_prepare_matrix <- function(se, assay_name, na_action = "complete",
                                 fill_value = -1, verbose = TRUE) {
  mat <- SummarizedExperiment::assay(se, assay_name)

  n_total <- nrow(mat)
  n_na_rows <- sum(apply(mat, 1, function(x) any(is.na(x))))

  if (n_na_rows == 0) {
    if (verbose) message("  Matrix has no missing values (", n_total, " proteins).")
    return(mat)
  }

  na_action <- match.arg(na_action, c("complete", "fill", "none"))

  if (na_action == "complete") {
    complete_rows <- complete.cases(mat)
    mat <- mat[complete_rows, , drop = FALSE]
    if (verbose)
      message("  na_action='complete': kept ", nrow(mat), " / ", n_total,
              " proteins (removed ", n_na_rows, " with NAs).")
    if (nrow(mat) < 10)
      warning("Very few proteins remaining (", nrow(mat),
              "). Consider using na_action='fill' or an imputed assay.")
  } else if (na_action == "fill") {
    mat[is.na(mat)] <- fill_value
    if (verbose)
      message("  na_action='fill': replaced NAs with ", fill_value,
              " (", n_total, " proteins).")
  } else {
    if (any(is.na(mat)))
      stop("Matrix contains NAs and na_action='none'. ",
           "Use na_action='complete' or 'fill', or provide a complete matrix.")
  }

  mat
}

#' Check which interaction terms are safe (n_levels < n_obs)
#'
#' Returns only interactions whose combined factor levels are strictly less
#' than the number of observations. This prevents lme4 from failing on
#' saturated random effects (e.g., Condition:Patient in paired designs).
#'
#' @param annot_df data.frame with factor columns
#' @param factors Character vector of main effect names
#' @param verbose Logical
#' @return Character vector of safe interaction terms (e.g., "A:B")
#' @keywords internal
.pvca_safe_interactions <- function(annot_df, factors, verbose = TRUE) {
  n_obs <- nrow(annot_df)
  if (length(factors) < 2) return(character(0))

  pairs <- combn(factors, 2, simplify = FALSE)
  safe <- character(0)
  skipped <- character(0)

  for (p in pairs) {
    interaction_levels <- nlevels(interaction(annot_df[[p[1]]],
                                              annot_df[[p[2]]],
                                              drop = TRUE))
    label <- paste(p, collapse = ":")
    if (interaction_levels < n_obs) {
      safe <- c(safe, label)
    } else {
      skipped <- c(skipped, label)
    }
  }

  if (verbose && length(skipped) > 0)
    message("  Skipped interaction(s) with n_levels >= n_obs: ",
            paste(skipped, collapse = ", "))

  safe
}

#' Core PVCA algorithm using lme4
#'
#' 1. PCA on the expression matrix
#' 2. Retain PCs explaining >= pca_threshold cumulative variance
#' 3. For each retained PC, fit lme4 mixed model with factors as random effects
#' 4. Extract variance components, weight by PC variance proportion
#'
#' @param mat Numeric matrix (proteins x samples)
#' @param annot_df data.frame with factor columns (rows = samples)
#' @param factors Character vector of factor names
#' @param pca_threshold Numeric (0-1)
#' @param verbose Logical
#' @return data.frame with columns: label, weights
#' @keywords internal
.pvca_run <- function(mat, annot_df, factors, pca_threshold = 0.6,
                      verbose = TRUE) {

  # --- PCA ---
  pca_res <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_pct <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  cum_var <- cumsum(var_pct)
  n_pcs <- min(which(cum_var >= pca_threshold))
  # Ensure at least 1 PC
  n_pcs <- max(1, n_pcs)
  # Cap at available PCs
  n_pcs <- min(n_pcs, length(var_pct))

  if (verbose)
    message("  PCA: retaining ", n_pcs, " PC(s) (",
            round(cum_var[n_pcs] * 100, 1), "% cumulative variance)")

  pc_scores <- pca_res$x[, seq_len(n_pcs), drop = FALSE]
  pc_weights <- var_pct[seq_len(n_pcs)]
  # Renormalize weights to sum to 1
  pc_weights <- pc_weights / sum(pc_weights)

  # --- Determine safe interactions ---
  safe_interactions <- .pvca_safe_interactions(annot_df, factors, verbose)

  # --- Build random-effects terms ---
  all_terms <- c(factors, safe_interactions)
  re_formula_str <- paste0("(1|", all_terms, ")", collapse = " + ")
  # Full formula: y ~ 1 + (1|factor1) + ... + (1|factorA:factorB) + ...

  # --- Fit mixed model for each PC and extract variance components ---
  # Prepare data frame for lme4
  fit_df <- as.data.frame(annot_df[, factors, drop = FALSE])
  # Ensure all are factors
  for (f in factors) fit_df[[f]] <- as.factor(fit_df[[f]])

  # Initialize accumulator for weighted variance proportions
  varcomp_accum <- setNames(rep(0, length(all_terms) + 1),
                            c(all_terms, "resid"))
  # Peso efectivo de los PCs realmente usados (para renormalizar si alguno se
  # salta por fallo de lmer o varianza total nula; asi los pesos suman 1).
  used_weight <- 0

  for (k in seq_len(n_pcs)) {
    fit_df$y <- pc_scores[, k]
    formula_str <- paste0("y ~ 1 + ", re_formula_str)
    # tryCatch: un fallo de convergencia en un solo PC no debe abortar toda la
    # descomposicion de varianza; se salta ese PC y se sigue con el resto.
    fm <- tryCatch(
      lme4::lmer(as.formula(formula_str), data = fit_df,
                 REML = TRUE,
                 control = lme4::lmerControl(
                   check.nobs.vs.nlev  = "warning",
                   check.nobs.vs.nRE   = "warning",
                   check.nlev.gtr.1    = "warning"
                 )),
      error = function(e) {
        if (verbose)
          message("  PVCA: lmer fallo en PC", k, " (", conditionMessage(e),
                  "); se omite este PC.")
        NULL
      }
    )
    if (is.null(fm)) next

    # Extract variance components
    vc <- lme4::VarCorr(fm)
    vc_df <- as.data.frame(vc)
    # vc_df has columns: grp, var1, var2, vcov, sdcor
    # grp contains factor names and "Residual"

    total_var <- sum(vc_df$vcov)
    if (total_var <= 0) next

    for (row_i in seq_len(nrow(vc_df))) {
      grp <- vc_df$grp[row_i]
      prop <- vc_df$vcov[row_i] / total_var
      key <- if (grp == "Residual") "resid" else grp
      if (key %in% names(varcomp_accum)) {
        varcomp_accum[key] <- varcomp_accum[key] + prop * pc_weights[k]
      }
    }
    used_weight <- used_weight + pc_weights[k]
  }

  # Renormalizar sobre el peso efectivo de los PCs usados para que los pesos
  # sumen 1 aunque se hayan omitido PCs.
  if (used_weight > 0 && used_weight < 1) {
    varcomp_accum <- varcomp_accum / used_weight
  }

  # --- Build result data.frame ---
  data.frame(
    label   = names(varcomp_accum),
    weights = as.numeric(varcomp_accum),
    stringsAsFactors = FALSE
  )
}

#' Categorize PVCA labels as technical, biological, interaction, or residual
#' @keywords internal
.pvca_categorize <- function(pvca_df, technical_factors, biological_factors,
                             variance_threshold = 0.01) {
  # Build interaction label sets
  tech_interactions <- if (length(technical_factors) >= 2) {
    g <- expand.grid(technical_factors, technical_factors, stringsAsFactors = FALSE)
    unique(paste(g$Var1, g$Var2, sep = ":"))
  } else character(0)

  biol_interactions <- if (length(biological_factors) >= 2) {
    g <- expand.grid(biological_factors, biological_factors, stringsAsFactors = FALSE)
    unique(paste(g$Var1, g$Var2, sep = ":"))
  } else character(0)

  label_of_small <- sprintf("Below %1.0f%%", 100 * variance_threshold)
  all_tech <- c(technical_factors, tech_interactions)
  all_biol <- c(biological_factors, biol_interactions)

  pvca_df$category <- ifelse(
    pvca_df$label %in% all_tech, "technical",
    ifelse(pvca_df$label %in% all_biol, "biological",
      ifelse(pvca_df$label %in% c(label_of_small, "resid"), "residual",
        "biol:techn"
      )
    )
  )

  # Sort: by weight descending, then small and resid at the end
  pvca_df <- pvca_df[order(-pvca_df$weights), ]
  small_idx <- pvca_df$label == label_of_small
  resid_idx <- pvca_df$label == "resid"
  other_idx <- !(small_idx | resid_idx)
  pvca_df <- rbind(
    pvca_df[other_idx, ],
    pvca_df[small_idx, ],
    pvca_df[resid_idx, ]
  )
  rownames(pvca_df) <- NULL
  pvca_df
}

#' Export a ggplot2 plot to file
#' @keywords internal
.pvca_export_gg_plot <- function(gg, filepath, width = 10, height = 6,
                                 dpi = 150) {
  ggplot2::ggsave(filename = filepath, plot = gg,
                  width = width, height = height, dpi = dpi)
}

# =============================================================================
# SECTION 2: PUBLIC FUNCTIONS
# =============================================================================

#' Compute PVCA variance components
#'
#' Performs Principal Variance Component Analysis on a SummarizedExperiment
#' object. Decomposes total variance into contributions from each specified
#' factor, their pairwise interactions, and residual variance.
#'
#' Uses lme4 directly to fit mixed models, which allows automatic detection
#' and exclusion of interaction terms where the number of levels equals or
#' exceeds the number of observations (e.g., Condition:Patient in paired
#' designs with one observation per patient-condition combination).
#'
#' @param se SummarizedExperiment object.
#' @param assay_name Character. Name of the assay to use. If NULL, uses the
#'   second assay (typically the normalized one).
#' @param factors Character vector. Column names in colData(se) to include as
#'   variance components. All specified factors must exist in colData.
#' @param pca_threshold Numeric (0-1). Minimum cumulative proportion of
#'   variance explained by retained principal components (default 0.6).
#' @param variance_threshold Numeric (0-1). Factors explaining less than this
#'   proportion are grouped into a "Below X%" category (default 0.01).
#' @param na_action Character. How to handle NAs: "complete" (default) removes
#'   rows with any NA, "fill" replaces NAs with fill_value, "none" assumes no
#'   NAs (errors if present).
#' @param fill_value Numeric. Replacement for NAs when na_action = "fill"
#'   (default -1).
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return data.frame with columns:
#'   \describe{
#'     \item{label}{Factor name, interaction term, "Below X%", or "resid"}
#'     \item{weights}{Weighted average proportion of variance (0-1)}
#'   }
#'
#' @examples
#' \dontrun{
#' vc <- pvca_compute(se, assay_name = "cycloess",
#'                    factors = c("Condition", "Gender", "Patient", "Injection"))
#' }
#' @export
pvca_compute <- function(se,
                         assay_name          = NULL,
                         factors,
                         pca_threshold       = 0.6,
                         variance_threshold  = 0.01,
                         na_action           = "complete",
                         fill_value          = -1,
                         verbose             = TRUE) {

  # --- Check required packages ---
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for pvca_compute(). ",
         "Install with: install.packages('lme4')")

  # --- Resolve assay name ---
  all_assays <- SummarizedExperiment::assayNames(se)
  if (is.null(assay_name)) {
    assay_name <- if (length(all_assays) >= 2) all_assays[2] else all_assays[1]
  }
  if (!assay_name %in% all_assays)
    stop("Assay '", assay_name, "' not found. Available: ",
         paste(all_assays, collapse = ", "))

  if (verbose) message("PVCA: using assay '", assay_name, "'")

  # --- Validate factors ---
  .pvca_validate_factors(se, factors)

  # --- Prepare matrix (handle NAs) ---
  mat <- .pvca_prepare_matrix(se, assay_name, na_action, fill_value, verbose)

  # --- Prepare annotation data ---
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  annot_df <- cd[colnames(mat), factors, drop = FALSE]

  # --- Run PVCA via lme4 ---
  if (verbose) message("  Running PVCA (", length(factors),
                       " factors, pca_threshold=", pca_threshold, ") ...")
  pvca_df <- .pvca_run(mat, annot_df, factors, pca_threshold, verbose)

  # Aggregate small components
  label_of_small <- sprintf("Below %1.0f%%", 100 * variance_threshold)
  small_mask <- pvca_df$weights < variance_threshold
  if (sum(small_mask) > 1) {
    pvca_small <- sum(pvca_df$weights[small_mask])
    pvca_df <- pvca_df[!small_mask, , drop = FALSE]
    pvca_df <- rbind(pvca_df, data.frame(label = label_of_small,
                                         weights = pvca_small))
  }

  if (verbose) {
    message("  PVCA complete. Variance components:")
    top <- pvca_df[order(-pvca_df$weights), ]
    for (i in seq_len(nrow(top))) {
      message(sprintf("    %-25s %5.1f%%", top$label[i], top$weights[i] * 100))
    }
  }

  pvca_df
}


#' Plot PVCA variance components
#'
#' Creates a bar plot showing the weighted average proportion of variance
#' attributed to each factor, colored by category (technical, biological,
#' interaction, residual).
#'
#' @param pvca_res data.frame. Output from pvca_compute() with columns
#'   'label' and 'weights'. Optionally pre-categorized with 'category'.
#' @param technical_factors Character vector. Factor names classified as
#'   technical. Required if pvca_res lacks a 'category' column.
#' @param biological_factors Character vector. Factor names classified as
#'   biological. Required if pvca_res lacks a 'category' column.
#' @param variance_threshold Numeric. Used for categorization if pvca_res
#'   lacks 'category' (default 0.01).
#' @param colors Named character vector of 4 colors for categories:
#'   "residual", "biological", "biol:techn", "technical". If NULL, defaults
#'   are used.
#' @param title Character. Plot title (default NULL).
#' @param base_size Numeric. Base font size for theme_classic (default 15).
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' vc <- pvca_compute(se, factors = c("Condition", "Gender", "Injection"))
#' pvca_plot(vc, technical_factors = "Injection",
#'           biological_factors = c("Condition", "Gender"))
#' }
#' @export
pvca_plot <- function(pvca_res,
                      technical_factors   = NULL,
                      biological_factors  = NULL,
                      variance_threshold  = 0.01,
                      colors              = NULL,
                      title               = NULL,
                      base_size           = 15) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for pvca_plot().")

  # --- Categorize if needed ---
  if (!"category" %in% colnames(pvca_res)) {
    if (is.null(technical_factors) && is.null(biological_factors))
      stop("pvca_res has no 'category' column. Provide technical_factors ",
           "and/or biological_factors for categorization.")
    tech <- technical_factors %||% character(0)
    biol <- biological_factors %||% character(0)
    pvca_res <- .pvca_categorize(pvca_res, tech, biol, variance_threshold)
  }

  # --- Factor ordering (preserve current row order) ---
  pvca_res$label <- factor(pvca_res$label, levels = pvca_res$label)

  # --- Colors (consistent with benchmarking palette) ---
  if (is.null(colors)) {
    colors <- c(
      "biological" = "#1a9850",
      "technical"  = "#E63946",
      "biol:techn" = "#F4D35E",
      "residual"   = "#CCCCCC"
    )
  }

  # --- Percentage labels for bars ---
  pvca_res$pct_label <- sprintf("%.1f%%", pvca_res$weights * 100)

  # --- Build plot ---
  gg <- ggplot2::ggplot(pvca_res,
                        ggplot2::aes(x = label, y = weights, fill = category)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = pct_label),
      vjust = -0.5, size = 3, color = "#495057"
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(x * 100), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.10))
    ) +
    ggplot2::labs(y = "Weighted Average Proportion of Variance", x = NULL,
                  fill = "Category") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 15, color = "#1D3557"),
      axis.text.x   = ggplot2::element_text(
        size = 11, angle = 45, hjust = 1, color = "#495057"),
      axis.text.y   = ggplot2::element_text(size = 11, color = "#495057"),
      axis.title.y  = ggplot2::element_text(size = 12, color = "#1D3557"),
      legend.position = "bottom",
      legend.title  = ggplot2::element_text(face = "bold"),
      legend.text   = ggplot2::element_text(size = 10),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank()
    )

  if (!is.null(title))
    gg <- gg + ggplot2::ggtitle(title)

  gg
}


#' PVCA Analysis — Orchestrator
#'
#' Runs the full PVCA workflow: compute variance components, generate plot,
#' and optionally export results.
#'
#' @param se SummarizedExperiment object.
#' @param assay_name Character. Assay to analyze. NULL = second assay
#'   (normalized, pre-imputation).
#' @param technical_factors Character vector. Factors representing technical
#'   variation (e.g., "Injection", "Digestion").
#' @param biological_factors Character vector. Factors representing biological
#'   variation (e.g., "Condition", "Gender", "Age").
#' @param pca_covariates Character vector or NULL. If provided, generates PCA
#'   plots colored by each covariate (grid + individual plots). These are
#'   included in the returned list and exported as PNGs. Default NULL (skip).
#' @param pca_threshold Numeric (0-1). Cumulative variance threshold for PCA
#'   (default 0.6).
#' @param variance_threshold Numeric (0-1). Minimum weight for individual
#'   display (default 0.01).
#' @param na_action Character: "complete", "fill", or "none" (default
#'   "complete").
#' @param fill_value Numeric. NA replacement when na_action = "fill"
#'   (default -1).
#' @param colors Named character vector for category colors (default NULL =
#'   built-in palette).
#' @param verbose Logical (default TRUE).
#' @param output_dir Character. Directory for exports. NULL = no export.
#' @param export_plots Logical (default TRUE).
#' @param export_tables Logical (default TRUE).
#' @param plot_width Numeric in inches (default 10).
#' @param plot_height Numeric in inches (default 6).
#' @param plot_dpi Numeric (default 150).
#'
#' @return Named list:
#'   \describe{
#'     \item{variance_components}{data.frame with label, weights, category}
#'     \item{plot}{ggplot2 object}
#'     \item{n_proteins}{Number of proteins used in the analysis}
#'     \item{assay_name}{Assay analyzed}
#'     \item{parameters}{List of parameters used}
#'   }
#'
#' @examples
#' \dontrun{
#' source("R/Batch_Correction.R")
#' pvca_res <- pvca_analysis(
#'   se = result$se_proc,
#'   assay_name = "cycloess",
#'   technical_factors = c("Injection", "Digestion"),
#'   biological_factors = c("Condition", "Gender", "Age", "Obesity",
#'                          "IMC", "Patient"),
#'   output_dir = "./results/my_analysis/"
#' )
#' pvca_res$variance_components
#' pvca_res$plot
#' }
#' @export
pvca_analysis <- function(se,
                          assay_name          = NULL,
                          technical_factors   = character(0),
                          biological_factors  = character(0),
                          pca_covariates      = NULL,
                          pca_threshold       = 0.6,
                          variance_threshold  = 0.01,
                          na_action           = "complete",
                          fill_value          = -1,
                          colors              = NULL,
                          verbose             = TRUE,
                          output_dir          = NULL,
                          export_plots        = TRUE,
                          export_tables       = TRUE,
                          plot_width          = 10,
                          plot_height         = 6,
                          plot_dpi            = 150) {

  # --- Required packages ---
  for (pkg in c("ggplot2", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for pvca_analysis().")
  }

  factors <- c(technical_factors, biological_factors)
  if (length(factors) == 0)
    stop("At least one factor must be specified via technical_factors ",
         "and/or biological_factors.")

  if (verbose) message("=== PVCA Analysis ===")

  # --- STEP 1: Compute variance components ---
  pvca_df <- pvca_compute(
    se                 = se,
    assay_name         = assay_name,
    factors            = factors,
    pca_threshold      = pca_threshold,
    variance_threshold = variance_threshold,
    na_action          = na_action,
    fill_value         = fill_value,
    verbose            = verbose
  )

  # Resolve actual assay_name used (for return)
  all_assays <- SummarizedExperiment::assayNames(se)
  used_assay <- assay_name %||%
    (if (length(all_assays) >= 2) all_assays[2] else all_assays[1])

  # Count proteins used
  mat <- .pvca_prepare_matrix(se, used_assay, na_action, fill_value,
                              verbose = FALSE)
  n_proteins <- nrow(mat)

  # --- STEP 2: Categorize ---
  pvca_df <- .pvca_categorize(pvca_df, technical_factors, biological_factors,
                              variance_threshold)

  # --- STEP 3: Plot ---
  if (verbose) message("  Generating PVCA plot ...")
  gg <- tryCatch(
    pvca_plot(pvca_df, colors = colors,
              title = paste0("PVCA \u2014 ", used_assay, " (", n_proteins,
                             " proteins)"),
              base_size = 15),
    error = function(e) {
      warning("pvca_plot() failed: ", conditionMessage(e))
      NULL
    }
  )

  # --- STEP 3b: PCA colored by covariates (optional) ---
  pca_cov_result <- NULL
  if (!is.null(pca_covariates) && length(pca_covariates) > 0) {
    if (verbose) message("  Generating PCA covariate plots ...")
    pca_cov_result <- tryCatch(
      pca_covariates_plot(
        se         = se,
        assay_name = used_assay,
        covariates = pca_covariates,
        na_action  = na_action,
        fill_value = fill_value,
        verbose    = verbose
      ),
      error = function(e) {
        warning("pca_covariates_plot() failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  # --- STEP 4: Export ---
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    if (export_tables) {
      tsv_path <- file.path(output_dir, "pvca_variance_components.tsv")
      if (requireNamespace("readr", quietly = TRUE)) {
        readr::write_tsv(pvca_df, tsv_path)
      } else {
        write.table(pvca_df, tsv_path, sep = "\t", row.names = FALSE,
                    quote = FALSE)
      }
      if (verbose) message("  Exported: ", tsv_path)
    }

    if (export_plots && !is.null(gg)) {
      png_path <- file.path(output_dir, "pvca_plot.png")
      .pvca_export_gg_plot(gg, png_path, plot_width, plot_height, plot_dpi)
      if (verbose) message("  Exported: ", png_path)
    }

    # Export PCA covariate plots
    if (export_plots && !is.null(pca_cov_result)) {
      grid_path <- file.path(output_dir, "pca_covariates_grid.png")
      .pvca_export_gg_plot(pca_cov_result$grid, grid_path,
                           plot_width, plot_height + 2, plot_dpi)
      if (verbose) message("  Exported: ", grid_path)

      for (cov_name in names(pca_cov_result$plots)) {
        cov_path <- file.path(output_dir,
                              paste0("pca_covariate_", cov_name, ".png"))
        .pvca_export_gg_plot(pca_cov_result$plots[[cov_name]], cov_path,
                             plot_width, plot_height, plot_dpi)
      }
      if (verbose) message("  Exported: ", length(pca_cov_result$plots),
                           " individual PCA covariate plots")
    }
  }

  if (verbose) message("=== PVCA Analysis complete ===")

  list(
    variance_components = pvca_df,
    plot                = gg,
    pca_covariates      = pca_cov_result,
    n_proteins          = n_proteins,
    assay_name          = used_assay,
    parameters          = list(
      technical_factors  = technical_factors,
      biological_factors = biological_factors,
      pca_covariates     = pca_covariates,
      pca_threshold      = pca_threshold,
      variance_threshold = variance_threshold,
      na_action          = na_action,
      fill_value         = fill_value
    )
  )
}


# =============================================================================
# SECTION 3: BATCH CORRECTION — INTERNAL HELPERS
# =============================================================================

#' Check that BERT is available
#' @keywords internal
.bc_check_bert <- function() {
  if (!requireNamespace("BERT", quietly = TRUE))
    stop("Package 'BERT' is required for batch correction.\n",
         "  Install with: BiocManager::install('BERT')")
  invisible(TRUE)
}

#' Validate batch column in colData
#' @param se SummarizedExperiment
#' @param batch_column Character. Column name in colData
#' @keywords internal
.bc_validate_batch <- function(se, batch_column) {
  cd_cols <- colnames(SummarizedExperiment::colData(se))
  if (!batch_column %in% cd_cols)
    stop("Batch column '", batch_column, "' not found in colData(se).\n",
         "  Available columns: ", paste(cd_cols, collapse = ", "), "\n",
         "  Ensure covariate_df with a '", batch_column,
         "' column is passed to process_proteomics().")

  batch_vals <- SummarizedExperiment::colData(se)[[batch_column]]
  n_na <- sum(is.na(batch_vals))
  if (n_na > 0)
    stop("Batch column '", batch_column, "' has ", n_na,
         " sample(s) with NA batch. Assign a batch to every sample (or remove ",
         "those samples) before batch correction — BERT/ComBat cannot handle ",
         "NA batch labels.")

  n_batch <- length(unique(batch_vals))
  if (n_batch < 2)
    stop("Batch column '", batch_column, "' has ", n_batch,
         " unique value(s). Batch correction requires at least 2 batches.")

  invisible(TRUE)
}

#' Identify features ComBat cannot fit
#'
#' Flags rows that, in any batch, have fewer than 2 finite observations or zero
#' within-batch variance. ComBat's parametric empirical-Bayes estimation divides
#' by the within-batch variance, so such features yield NaN and abort the whole
#' correction (`while (change > conv)` receives NA). These features cannot be
#' batch-corrected anyway (constant within a batch) and should pass through
#' unadjusted.
#'
#' @param mat Numeric matrix (features x samples).
#' @param batch_vec Batch assignment per column.
#' @return Logical vector (length nrow(mat)); TRUE = unfittable by ComBat.
#' @keywords internal
.bc_combat_unfittable <- function(mat, batch_vec) {
  bad <- logical(nrow(mat))
  for (b in unique(batch_vec[!is.na(batch_vec)])) {
    sub   <- mat[, batch_vec == b, drop = FALSE]
    n_obs <- rowSums(is.finite(sub))
    rvar  <- apply(sub, 1, function(z) {
      z <- z[is.finite(z)]
      if (length(z) < 2L) NA_real_ else stats::var(z)
    })
    bad <- bad | n_obs < 2L | is.na(rvar) | rvar == 0
  }
  bad
}

#' Run BERT batch correction
#'
#' Core wrapper around BERT::BERT() with error handling.
#' BERT uses hierarchical tree decomposition: pairwise batch corrections
#' are organized in a binary tree, enabling robust handling of missing values
#' and multiple batches.
#'
#' @param mat Numeric matrix (proteins x samples)
#' @param batch_vec Factor or character vector of batch assignments per sample
#' @param method Character: "ComBat" (default), "limma", or "ref"
#' @param combatmode Integer 1-4 (only used when method = "ComBat")
#' @param covariates data.frame with columns Cov_1, Cov_2, ... or NULL
#' @param qualitycontrol Logical: compute ASW quality metrics
#' @return Numeric matrix (proteins x samples), possibly fewer rows
#' @keywords internal
#' @noRd
.bc_run_bert <- function(mat, batch_vec,
                         method         = "ComBat",
                         combatmode     = 1,
                         covariates     = NULL,
                         qualitycontrol = FALSE) {

  orig_colnames <- colnames(mat)
  n_samples  <- ncol(mat)
  n_features <- nrow(mat)

  # Transpose: proteins x samples → samples x features (BERT format)
  bert_input <- as.data.frame(t(mat))

  # Add Batch column (integer required by BERT)
  bert_input$Batch <- as.integer(as.factor(batch_vec))

  # Add covariate columns if provided (Cov_1, Cov_2, ...)
  if (!is.null(covariates)) {
    # Realinear por nombre de muestra antes del cbind (defensivo: no depender
    # de que el orden de filas de covariates coincida con colnames(mat)).
    if (!is.null(rownames(covariates)) &&
        all(orig_colnames %in% rownames(covariates))) {
      covariates <- covariates[orig_colnames, , drop = FALSE]
    }
    bert_input <- cbind(bert_input, covariates)
  }

  # Call BERT
  result <- suppressWarnings(
    BERT::BERT(bert_input,
               method         = method,
               combatmode     = combatmode,
               qualitycontrol = qualitycontrol)
  )

  # Extract corrected matrix (remove Batch, Cov_* y nombres reservados de BERT)
  meta_cols <- c("Batch", "Label", "Sample", "Reference", "Cov",
                 grep("^Cov_", colnames(result), value = TRUE))
  corrected <- as.matrix(result[, !colnames(result) %in% meta_cols,
                                drop = FALSE])

  # Transpose back: samples x features → proteins x samples
  corrected <- t(corrected)

  # Verify dimensions
  if (ncol(corrected) != n_samples) {
    stop("BERT returned ", ncol(corrected), " of ", n_samples,
         " samples. Check batch assignments and data structure.")
  }

  # Reorder to match input
  corrected <- corrected[, orig_colnames, drop = FALSE]

  corrected
}


# =============================================================================
# SECTION 4: BATCH CORRECTION — PUBLIC FUNCTION
# =============================================================================

#' Batch Correction with BERT
#'
#' Applies BERT batch effect correction to a normalized assay in a
#' SummarizedExperiment. BERT uses hierarchical tree decomposition to handle
#' missing values, making it suitable for pre-imputation batch correction on
#' protein-level proteomics data.
#'
#' @param se SummarizedExperiment with a normalized assay.
#' @param assay_name Character. Name of the input assay to correct
#'   (e.g., "cycloess").
#' @param batch_column Character. Column in colData(se) containing batch
#'   assignments (default "Batch").
#' @param corrected_assay_name Character. Name for the new corrected assay
#'   added to the SE (default "BERT").
#' @param algorithm Character: "ComBat" (default), "limma", or "ref".
#' @param ComBat_mode Integer 1-4 controlling ComBat behavior:
#'   1 = parametric + mean+variance, 2 = parametric + mean-only,
#'   3 = non-parametric + mean+variance, 4 = non-parametric + mean-only.
#' @param covariates Character vector of column names from colData(se) to
#'   use as categorical covariates for batch correction (default NULL).
#'   These are mapped to BERT's Cov_1, Cov_2, ... format internally.
#'   IMPORTANTE: para ComBat/limma, incluye aqui la variable biologica de
#'   interes (p.ej. la condicion) para PRESERVARLA; de lo contrario ComBat
#'   elimina toda la varianza del batch y puede borrar senal biologica si
#'   condicion y batch estan confundidos. Deben ser covariables categoricas.
#' @param qualitycontrol Logical. Compute ASW (Average Silhouette Width)
#'   quality metrics for raw vs corrected data (default FALSE).
#' @param verbose Logical (default TRUE).
#'
#' @return SummarizedExperiment with the new corrected assay added.
#'   If BERT returns fewer features than the input, the SE is subsetted
#'   to match and a warning is issued.
#'
#' @examples
#' \dontrun{
#' source("R/Batch_Correction.R")
#' se_corrected <- batch_correct_proteomics(
#'   se         = result$se_proc,
#'   assay_name = "cycloess",
#'   batch_column = "Batch"
#' )
#' SummarizedExperiment::assayNames(se_corrected)
#' # [1] "raw" "log2" "cycloess" "BERT"
#'
#' # With covariates
#' se_corrected <- batch_correct_proteomics(
#'   se         = result$se_proc,
#'   assay_name = "cycloess",
#'   batch_column = "Batch",
#'   covariates = c("Gender", "Age")
#' )
#' }
#' @export
batch_correct_proteomics <- function(
    se,
    assay_name,
    batch_column           = "Batch",
    corrected_assay_name   = "BERT",
    algorithm              = "ComBat",
    ComBat_mode            = 1,
    covariates             = NULL,
    qualitycontrol         = FALSE,
    verbose                = TRUE
) {

  # --- Check dependencies ---
  .bc_check_bert()

  # --- Validate algorithm ---
  algorithm <- match.arg(algorithm, c("ComBat", "limma", "ref"))

  # --- Validate assay ---
  all_assays <- SummarizedExperiment::assayNames(se)
  if (!assay_name %in% all_assays)
    stop("Assay '", assay_name, "' not found in SE. Available: ",
         paste(all_assays, collapse = ", "))

  # --- Validate batch column ---
  .bc_validate_batch(se, batch_column)

  batch_vals <- SummarizedExperiment::colData(se)[[batch_column]]
  n_batches <- length(unique(batch_vals[!is.na(batch_vals)]))
  if (verbose) {
    cat("\n=== BATCH CORRECTION (BERT) ===\n")
    cat("- Input assay:", assay_name, "\n")
    cat("- Batch column:", batch_column,
        "(", n_batches, "batches )\n")
    cat("- Algorithm:", algorithm)
    if (algorithm == "ComBat") cat(" (mode", ComBat_mode, ")")
    cat("\n")
    if (!is.null(covariates))
      cat("- Covariates:", paste(covariates, collapse = ", "), "\n")
  }

  # --- Extract matrix ---
  mat <- SummarizedExperiment::assay(se, assay_name)
  n_features_in <- nrow(mat)
  na_count_in  <- sum(is.na(mat))
  na_total_in  <- length(mat)
  na_pct_in    <- round(100 * na_count_in / na_total_in, 2)
  if (verbose) {
    cat("- Input features:", n_features_in,
        "| Samples:", ncol(mat), "\n")
    cat("- Input NAs:", na_count_in, "/", na_total_in,
        "(", na_pct_in, "% )\n")
  }

  # --- Build covariates data.frame in BERT format (Cov_1, Cov_2, ...) ---
  cov_df <- NULL
  if (!is.null(covariates)) {
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    missing_covs <- setdiff(covariates, colnames(cd))
    if (length(missing_covs) > 0)
      stop("Covariate column(s) not found in colData: ",
           paste(missing_covs, collapse = ", "))
    cov_df <- cd[, covariates, drop = FALSE]
    colnames(cov_df) <- paste0("Cov_", seq_along(covariates))
    rownames(cov_df) <- colnames(mat)  # para el realineado por nombre en .bc_run_bert
  } else if (algorithm %in% c("ComBat", "limma")) {
    warning("batch_correct_proteomics: 'covariates = NULL' con algorithm='",
            algorithm, "'. ComBat/limma eliminan TODA la varianza asociada al ",
            "batch; si la condicion biologica esta (parcialmente) confundida con ",
            "el batch, se perdera senal biologica real. Se recomienda pasar la ",
            "variable de condicion en 'covariates' (batch_covariates) para ",
            "preservarla.", call. = FALSE)
  }

  # --- Apart features ComBat cannot fit (zero within-batch variance) ---
  # ComBat (parametric EB) divides by within-batch variance; features that are
  # constant within a batch produce NaN and abort the whole correction. Such
  # features cannot be batch-corrected anyway, so set them aside and pass them
  # through unadjusted.
  combat_bad <- logical(nrow(mat))
  if (algorithm == "ComBat") combat_bad <- .bc_combat_unfittable(mat, batch_vals)
  if (any(combat_bad) && verbose)
    cat("- Features no corregibles (varianza intra-batch nula):", sum(combat_bad),
        "-> se mantienen sin ajustar\n")

  # --- Run BERT ---
  if (verbose) cat("- Running BERT ...\n")
  corrected_fit <- .bc_run_bert(
    mat            = mat[!combat_bad, , drop = FALSE],
    batch_vec      = batch_vals,
    method         = algorithm,
    combatmode     = ComBat_mode,
    covariates     = cov_df,
    qualitycontrol = qualitycontrol
  )

  if (any(combat_bad)) {
    passthrough <- mat[combat_bad, colnames(corrected_fit), drop = FALSE]
    corrected_mat <- rbind(corrected_fit, passthrough)
    corrected_mat <- corrected_mat[intersect(rownames(mat), rownames(corrected_mat)),
                                   , drop = FALSE]
  } else {
    corrected_mat <- corrected_fit
  }

  # --- Align output to SE ---
  n_features_out <- nrow(corrected_mat)

  # Handle potential feature loss
  if (n_features_out < n_features_in) {
    n_dropped <- n_features_in - n_features_out
    warning("BERT returned ", n_features_out, " features (", n_dropped,
            " dropped due to batch-specific missingness).",
            "\n  The SE will be subsetted to match.")

    # Subset SE to keep only features returned by BERT
    keep_features <- rownames(corrected_mat)
    se <- se[keep_features, ]
    if (verbose) cat("- Features after correction:", n_features_out,
                     " (", n_dropped, " dropped)\n")
  } else {
    if (verbose) cat("- Features after correction:", n_features_out,
                     " (none dropped)\n")
  }

  # --- Verify NA pattern ---
  na_count_out <- sum(is.na(corrected_mat))
  na_total_out <- length(corrected_mat)
  na_pct_out   <- round(100 * na_count_out / na_total_out, 2)
  if (verbose) cat("- Output NAs:", na_count_out, "/", na_total_out,
                   "(", na_pct_out, "% )\n")

  if (na_count_out != na_count_in && n_features_out == n_features_in) {
    warning("NA count changed after BERT: ", na_count_in, " -> ", na_count_out,
            ". Expected identical NA pattern when no features are dropped.")
  }

  # --- Add corrected assay to SE ---
  # Ensure exact row/column order matches SE (BERT may reorder features)
  corrected_mat <- corrected_mat[rownames(se), colnames(se), drop = FALSE]
  SummarizedExperiment::assay(se, corrected_assay_name) <- corrected_mat

  if (verbose) {
    cat("- New assay added: '", corrected_assay_name, "'\n", sep = "")
    cat("- Assays in SE:", paste(SummarizedExperiment::assayNames(se),
                                 collapse = ", "), "\n")
    cat("=== BATCH CORRECTION COMPLETADA ===\n")
  }

  se
}


# =============================================================================
# SECTION 5: PCA COLORED BY COVARIATES
# =============================================================================

#' Compute PCA scores from a SE assay
#'
#' @param se SummarizedExperiment
#' @param assay_name Character. Assay to use.
#' @param na_action Character: "complete" (remove rows with NAs) or "fill".
#' @param fill_value Numeric. Replacement for NAs when na_action="fill".
#' @param center Logical. Center variables before PCA (default TRUE).
#' @param scale. Logical. Scale variables before PCA (default TRUE).
#' @param verbose Logical.
#' @return List with: scores (data.frame PC1, PC2, Sample + colData),
#'   pct_var (numeric vector of % variance per PC), n_proteins (integer).
#' @keywords internal
.pca_compute_scores <- function(se, assay_name, na_action = "complete",
                                fill_value = -1, center = TRUE,
                                scale. = TRUE, verbose = TRUE) {

  mat <- .pvca_prepare_matrix(se, assay_name, na_action, fill_value, verbose)

  if (nrow(mat) < 2)
    stop("Insufficient complete rows for PCA (", nrow(mat), " rows).")

  pca_res <- stats::prcomp(t(mat), center = center, scale. = scale.)
  pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  scores <- data.frame(
    PC1    = pca_res$x[, 1],
    PC2    = pca_res$x[, 2],
    Sample = rownames(pca_res$x),
    stringsAsFactors = FALSE
  )

  # Merge all colData
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cd$Sample <- colnames(se)
  scores <- merge(scores, cd, by = "Sample")

  list(scores = scores, pct_var = pct_var, n_proteins = nrow(mat))
}


#' PCA Plots Colored by Covariates
#'
#' Generates PCA scatter plots (PC1 vs PC2) from a SummarizedExperiment assay,
#' coloring points by different covariates. Useful for visually identifying
#' batch effects or confounding.
#'
#' Categorical covariates use a discrete color scale; numeric covariates
#' use a continuous viridis gradient.
#'
#' @param se SummarizedExperiment object.
#' @param assay_name Character. Assay to analyze. NULL = second assay
#'   (normalized, pre-imputation).
#' @param covariates Character vector. Column names in colData(se) to use
#'   for coloring (e.g., c("Injection", "Digestion", "Condition")).
#' @param na_action Character: "complete" (default) or "fill".
#' @param fill_value Numeric. NA replacement when na_action="fill" (default -1).
#' @param center Logical. Center before PCA (default TRUE).
#' @param scale. Logical. Scale before PCA (default TRUE).
#' @param de_results Data.frame with differential expression results (e.g.,
#'   result$DEPs_results). Must contain columns "Protein.IDs", "adj.P.Val", and
#'   "Comparison". When provided, PCA is computed using only significant proteins.
#'   Default NULL (all proteins).
#' @param comparison Character. Which comparison to filter by. Use "any"
#'   (default) to keep proteins significant in any comparison, or specify a
#'   comparison name (e.g., "Post_vs_Pre") to filter to that one only.
#'   Ignored if de_results is NULL.
#' @param alpha Numeric. Significance threshold for adj.P.Val when filtering
#'   by de_results (default 0.05). Ignored if de_results is NULL.
#' @param filter_samples Logical. If TRUE and a specific comparison is given
#'   (not "any"), subset samples to only those belonging to the comparison
#'   conditions. Mirrors pca_highchart_list(filter_samples_to_comparison).
#'   Default FALSE (use all samples).
#' @param condition_column Character. Column name in colData(se) containing
#'   condition labels, used when filter_samples=TRUE. Default "Condition".
#' @param point_size Numeric. Size of scatter points (default 3).
#' @param verbose Logical (default TRUE).
#'
#' @return Named list:
#'   \describe{
#'     \item{grid}{ggplot2 facet_wrap with one panel per covariate}
#'     \item{plots}{Named list of individual ggplot2 objects per covariate}
#'     \item{pca_summary}{data.frame with PC1/PC2 variance explained}
#'     \item{assay_name}{Assay used}
#'     \item{n_proteins}{Number of proteins used}
#'   }
#'
#' @examples
#' \dontrun{
#' source("R/Batch_Correction.R")
#' pca_cov <- pca_covariates_plot(
#'   se = result$se_proc,
#'   assay_name = "BERT",
#'   covariates = c("Injection", "Digestion", "Condition", "Gender")
#' )
#' pca_cov$grid              # faceted grid
#' pca_cov$plots$Injection   # individual plot
#' }
#' @export
pca_covariates_plot <- function(
    se,
    assay_name    = NULL,
    covariates,
    na_action     = "complete",
    fill_value    = -1,
    center        = TRUE,
    scale.        = TRUE,
    de_results       = NULL,
    comparison       = "any",
    alpha            = 0.05,
    filter_samples   = FALSE,
    condition_column = "Condition",
    point_size       = 3,
    verbose          = TRUE
) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for pca_covariates_plot().")

  # --- Resolve assay name ---
  all_assays <- SummarizedExperiment::assayNames(se)
  if (is.null(assay_name)) {
    assay_name <- if (length(all_assays) >= 2) all_assays[2] else all_assays[1]
  }
  if (!assay_name %in% all_assays)
    stop("Assay '", assay_name, "' not found. Available: ",
         paste(all_assays, collapse = ", "))

  # --- Validate covariates ---
  .pvca_validate_factors(se, covariates)

  # --- Filter to differential proteins if de_results provided ---
  protein_label <- "proteins"
  if (!is.null(de_results)) {
    required_cols <- c("Protein.IDs", "adj.P.Val", "Comparison")
    missing_cols <- setdiff(required_cols, colnames(de_results))
    if (length(missing_cols) > 0)
      stop("de_results must contain columns: ",
           paste(missing_cols, collapse = ", "))

    # Filter by comparison
    if (!identical(comparison, "any")) {
      available_comps <- unique(de_results$Comparison)
      if (!comparison %in% available_comps)
        stop("Comparison '", comparison, "' not found. Available: ",
             paste(available_comps, collapse = ", "))
      de_sub <- de_results[de_results$Comparison == comparison, ]
      comp_label <- comparison
    } else {
      de_sub <- de_results
      comp_label <- "any comparison"
    }

    sig_ids <- unique(de_sub$Protein.IDs[de_sub$adj.P.Val < alpha])
    keep <- rownames(se) %in% sig_ids
    if (sum(keep) < 3)
      stop("Only ", sum(keep), " differential proteins found (adj.P.Val < ",
           alpha, ", comparison: ", comp_label, "). Need at least 3 for PCA.")
    se <- se[keep, ]
    protein_label <- paste0("DEPs (", comp_label, ")")
    if (verbose) message("PCA covariates: filtered to ", sum(keep),
                         " DEPs (adj.P.Val < ", alpha,
                         ", comparison: ", comp_label, ")")
  }

  # --- Filter samples to comparison conditions ---
  if (isTRUE(filter_samples) && !identical(comparison, "any") &&
      !is.null(de_results)) {
    conds <- unique(trimws(unlist(strsplit(comparison, "[-|:_]"))))
    conds <- conds[nchar(conds) > 0 & conds != "vs"]

    if (!condition_column %in% colnames(SummarizedExperiment::colData(se)))
      stop("Column '", condition_column,
           "' not found in colData for sample filtering.")

    cd <- SummarizedExperiment::colData(se)
    keep_samples <- colnames(se)[cd[[condition_column]] %in% conds]

    if (length(keep_samples) < 3)
      stop("Only ", length(keep_samples),
           " samples after filtering by comparison '", comparison,
           "'. Need at least 3.")

    se <- se[, keep_samples]
    if (verbose) message("PCA covariates: filtered to ", length(keep_samples),
                         " samples from conditions: ",
                         paste(conds, collapse = ", "))
  }

  if (verbose) message("PCA covariates: using assay '", assay_name, "'")

  # --- Compute PCA ---
  pca_data <- .pca_compute_scores(
    se = se, assay_name = assay_name, na_action = na_action,
    fill_value = fill_value, center = center, scale. = scale.,
    verbose = verbose
  )

  scores     <- pca_data$scores
  pct_var    <- pca_data$pct_var
  n_proteins <- pca_data$n_proteins

  x_lab <- paste0("PC1 (", pct_var[1], "%)")
  y_lab <- paste0("PC2 (", pct_var[2], "%)")

  # --- Common theme (matching PCA_Highcharts style) ---
  common_theme <- ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 16, color = "#1D3557"),
      plot.subtitle    = ggplot2::element_text(
        hjust = 0.5, size = 10, color = "#495057"),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major = ggplot2::element_line(
        color = "#F1F3F4", linetype = "dotted", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(size = 11, color = "#495057"),
      axis.title       = ggplot2::element_text(
        size = 13, face = "bold", color = "#212529"),
      axis.line        = ggplot2::element_line(color = "#DEE2E6", linewidth = 0.5),
      axis.ticks       = ggplot2::element_line(color = "#DEE2E6", linewidth = 0.3),
      strip.text       = ggplot2::element_text(
        face = "bold", size = 12, color = "#1D3557"),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold", color = "#495057"),
      legend.text      = ggplot2::element_text(size = 10, color = "#495057"),
      plot.margin      = ggplot2::margin(10, 15, 10, 15)
    )

  # Reference lines at x=0 and y=0 (dashed, subtle)
  ref_lines <- list(
    ggplot2::geom_hline(yintercept = 0, color = "#ADB5BD",
                        linetype = "dashed", linewidth = 0.4),
    ggplot2::geom_vline(xintercept = 0, color = "#ADB5BD",
                        linetype = "dashed", linewidth = 0.4)
  )

  # Default palette: hcl Dark 3 (same as Highcharts PCA)
  .pca_cov_palette <- function(n) {
    grDevices::hcl.colors(max(n, 3), "Dark 3")[seq_len(n)]
  }

  # --- Individual plots ---
  plot_list <- stats::setNames(vector("list", length(covariates)), covariates)

  for (cov in covariates) {
    vals <- scores[[cov]]
    is_numeric <- is.numeric(vals) && length(unique(vals)) > 6

    if (!is_numeric) {
      scores[[cov]] <- as.factor(scores[[cov]])
    }

    gg <- ggplot2::ggplot(scores,
                          ggplot2::aes(x = PC1, y = PC2, color = .data[[cov]])) +
      ref_lines +
      ggplot2::geom_point(size = point_size, alpha = 0.85) +
      ggplot2::labs(
        title    = paste0("PCA \u2014 ", assay_name),
        subtitle = paste0(n_proteins, " ", protein_label, " | colored by ", cov),
        x = x_lab, y = y_lab, color = cov
      ) +
      common_theme

    if (is_numeric) {
      gg <- gg + ggplot2::scale_color_viridis_c(option = "D")
    } else {
      n_levels <- length(unique(vals))
      gg <- gg +
        ggplot2::scale_color_manual(values = .pca_cov_palette(n_levels))
    }

    plot_list[[cov]] <- gg
  }

  # --- Grid (faceted) ---
  # Build long-format data for faceting
  long_rows <- vector("list", length(covariates))
  for (i in seq_along(covariates)) {
    cov <- covariates[i]
    vals <- scores[[cov]]
    is_numeric <- is.numeric(vals) && length(unique(vals)) > 6
    long_rows[[i]] <- data.frame(
      PC1       = scores$PC1,
      PC2       = scores$PC2,
      Sample    = scores$Sample,
      Covariate = cov,
      Value     = if (is_numeric) as.character(round(as.numeric(vals), 1))
                  else as.character(vals),
      stringsAsFactors = FALSE
    )
  }
  long_df <- do.call(rbind, long_rows)
  long_df$Covariate <- factor(long_df$Covariate, levels = covariates)

  grid_plot <- ggplot2::ggplot(long_df,
                               ggplot2::aes(x = PC1, y = PC2, color = Value)) +
    ref_lines +
    ggplot2::geom_point(size = point_size * 0.7, alpha = 0.85) +
    ggplot2::facet_wrap(~ Covariate, ncol = 2) +
    ggplot2::scale_color_manual(
      values = .pca_cov_palette(length(unique(long_df$Value)))) +
    ggplot2::labs(
      title    = paste0("PCA \u2014 ", assay_name,
                         " (", n_proteins, " ", protein_label, ")"),
      x = x_lab, y = y_lab, color = "Value"
    ) +
    common_theme +
    ggplot2::theme(legend.position = "none")

  if (verbose) message("  Generated ", length(covariates),
                       " individual PCA plots + 1 grid plot.")

  list(
    grid        = grid_plot,
    plots       = plot_list,
    pca_summary = data.frame(
      PC       = paste0("PC", seq_along(pct_var)),
      Variance = pct_var
    ),
    assay_name  = assay_name,
    n_proteins  = n_proteins
  )
}
