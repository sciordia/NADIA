# =============================================================================
# Batch Diagnostics — PVCA + Batch Correction (HarmonizR)
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
# SECTION 3-4: Batch Correction (HarmonizR)
#   Optional batch effect correction using HarmonizR (ComBat/limma with
#   matrix dissection for missing-value-tolerant correction).
#   Public functions:
#     batch_correct_proteomics()   : Correct batch effects, add assay to SE
#
# References:
#   PVCA: Li et al. (2009) Biostatistics 10(2):317-326
#         Cuklina et al. (2021) Molecular & Cellular Proteomics 20 (proBatch)
#   HarmonizR: Voß et al. (2022) Nature Communications 13:6171
#              Gregoricchio et al. (2024) DEprot
#
# Implementation notes:
#   - PVCA uses lme4 directly (not pvca package) to allow automatic detection
#     and exclusion of problematic interaction terms where n_levels >= n_obs
#   - HarmonizR handles missing values via matrix dissection, suitable for
#     pre-imputation batch correction on protein-level data
#
# Dependencies:
#   PVCA: SummarizedExperiment, ggplot2, lme4
#   Batch Correction: SummarizedExperiment, HarmonizR (Bioconductor)
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# --- Null coalescing operator ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

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

  for (k in seq_len(n_pcs)) {
    fit_df$y <- pc_scores[, k]
    formula_str <- paste0("y ~ 1 + ", re_formula_str)
    fm <- lme4::lmer(as.formula(formula_str), data = fit_df,
                     REML = TRUE,
                     control = lme4::lmerControl(
                       check.nobs.vs.nlev  = "warning",
                       check.nobs.vs.nRE   = "warning",
                       check.nlev.gtr.1    = "warning"
                     ))

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

  # --- Colors ---
  if (is.null(colors)) {
    colors <- c(
      "residual"   = "#999999",
      "biological" = "#56B4E9",
      "biol:techn" = "#E69F00",
      "technical"  = "#D55E00"
    )
  }

  # --- Build plot ---
  gg <- ggplot2::ggplot(pvca_res,
                        ggplot2::aes(x = label, y = weights, fill = category)) +
    ggplot2::geom_bar(stat = "identity", color = "black", width = 0.7) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(y = "Weighted average proportion of variance", x = NULL,
                  fill = "Category") +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "top"
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
#' source("R/Batch_Diagnostics.R")
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
  }

  if (verbose) message("=== PVCA Analysis complete ===")

  list(
    variance_components = pvca_df,
    plot                = gg,
    n_proteins          = n_proteins,
    assay_name          = used_assay,
    parameters          = list(
      technical_factors  = technical_factors,
      biological_factors = biological_factors,
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

#' Check that HarmonizR is available
#' @keywords internal
.bc_check_harmonizr <- function() {
  if (!requireNamespace("HarmonizR", quietly = TRUE))
    stop("Package 'HarmonizR' is required for batch correction.\n",
         "  Install with: BiocManager::install('HarmonizR')")
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
  n_batch <- length(unique(batch_vals[!is.na(batch_vals)]))
  if (n_batch < 2)
    stop("Batch column '", batch_column, "' has ", n_batch,
         " unique value(s). Batch correction requires at least 2 batches.")

  invisible(TRUE)
}

#' Build HarmonizR description data.frame from SE colData
#'
#' Creates the batch description format required by HarmonizR:
#' a data.frame with sample IDs and numeric batch assignments.
#'
#' @param se SummarizedExperiment
#' @param batch_column Character. Column name containing batch info
#' @return data.frame with columns: sample (character), batch (integer)
#' @keywords internal
.bc_build_description <- function(se, batch_column) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  data.frame(
    ID     = colnames(se),
    sample = seq_len(ncol(se)),
    batch  = as.integer(as.factor(cd[[batch_column]])),
    stringsAsFactors = FALSE
  )
}

#' Run HarmonizR batch correction
#'
#' Core wrapper around HarmonizR::harmonizR() with error handling.
#'
#' @param mat Numeric matrix (proteins x samples)
#' @param description data.frame from .bc_build_description()
#' @param algorithm Character: "ComBat" or "limma"
#' @param ComBat_mode Integer 1-4
#' @param sort_method Character: sorting strategy for matrix dissection
#' @param block Integer or NULL
#' @param cores Integer
#' @param ur Logical: unique combination removal
#' @return Numeric matrix (proteins x samples), possibly fewer rows
#' @keywords internal
.bc_run_harmonizr <- function(mat, description,
                              algorithm   = "ComBat",
                              ComBat_mode = 1,
                              sort_method = "sparsity_sort",
                              block       = NULL,
                              cores       = 1,
                              ur          = TRUE) {

  # Use safe column/row names (S1..Sn, F1..Fm) to avoid any name mangling
  # through HarmonizR's file I/O pipeline. Restore originals after.
  orig_colnames <- colnames(mat)
  orig_rownames <- rownames(mat)
  n_samples  <- ncol(mat)
  n_features <- nrow(mat)

  safe_col <- paste0("S", seq_len(n_samples))
  safe_row <- paste0("F", seq_len(n_features))

  colnames(mat) <- safe_col
  rownames(mat) <- safe_row

  # Update description to use safe sample names
  description$ID <- safe_col

  # Create temp directory for HarmonizR I/O
  tmp_dir <- tempfile("harmonizr_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Write data file (TSV: first column = feature IDs, rest = samples)
  data_file <- file.path(tmp_dir, "input_data.tsv")
  data_df <- data.frame(ID = safe_row, mat, check.names = FALSE)
  utils::write.table(data_df, data_file, sep = "\t", quote = FALSE,
                      row.names = FALSE)

  # Write description file (CSV: ID, sample, batch)
  desc_file <- file.path(tmp_dir, "description.csv")
  utils::write.csv(description, desc_file, row.names = FALSE)

  # Output file path (HarmonizR appends .tsv)
  output_base <- file.path(tmp_dir, "cured_data")

  # sort is only useful with block; HarmonizR 1.8.0 has a sorting bug
  use_sort <- if (!is.null(block)) sort_method else FALSE

  hr_args <- list(
    data_as_input        = data_file,
    description_as_input = desc_file,
    algorithm            = algorithm,
    ComBat_mode          = ComBat_mode,
    sort                 = use_sort,
    cores                = cores,
    ur                   = ur,
    output_file          = output_base
  )
  if (!is.null(block)) hr_args$block <- block

  do.call(HarmonizR::harmonizR, hr_args)

  # Read the output file written by HarmonizR
  output_file <- paste0(output_base, ".tsv")
  if (!file.exists(output_file))
    stop("HarmonizR did not produce output file: ", output_file)

  result_df <- utils::read.delim(output_file, sep = "\t", row.names = 1,
                                  check.names = FALSE)
  result <- as.matrix(result_df)

  # Restore original names via safe→original mapping
  col_map <- setNames(orig_colnames, safe_col)
  row_map <- setNames(orig_rownames, safe_row)

  colnames(result) <- col_map[colnames(result)]
  rownames(result) <- row_map[rownames(result)]

  # Check for sample loss (too many batches + missing data can cause this)
  if (ncol(result) < n_samples) {
    stop("HarmonizR returned only ", ncol(result), " of ", n_samples,
         " samples.\n",
         "  This happens when matrix dissection cannot cover all samples ",
         "(too many batches relative to sample size and missingness).\n",
         "  Try: (1) fewer batches (e.g., correct by a single factor), ",
         "(2) set block parameter, or (3) algorithm='limma'.")
  }

  # Reorder to match input
  result <- result[, orig_colnames, drop = FALSE]

  result
}


# =============================================================================
# SECTION 4: BATCH CORRECTION — PUBLIC FUNCTION
# =============================================================================

#' Batch Correction with HarmonizR
#'
#' Applies HarmonizR batch effect correction to a normalized assay in a
#' SummarizedExperiment. HarmonizR uses matrix dissection to handle missing
#' values, making it suitable for pre-imputation batch correction on
#' protein-level proteomics data.
#'
#' @param se SummarizedExperiment with a normalized assay.
#' @param assay_name Character. Name of the input assay to correct
#'   (e.g., "cycloess").
#' @param batch_column Character. Column in colData(se) containing batch
#'   assignments (default "Batch").
#' @param corrected_assay_name Character. Name for the new corrected assay
#'   added to the SE (default "HarmonizR").
#' @param algorithm Character: "ComBat" (default) or "limma".
#' @param ComBat_mode Integer 1-4 controlling ComBat behavior:
#'   1 = parametric + mean+variance, 2 = parametric + mean-only,
#'   3 = non-parametric + mean+variance, 4 = non-parametric + mean-only.
#' @param sort_method Character. Sorting strategy for matrix dissection:
#'   "sparsity_sort" (default), "seriation_sort", or "jaccard_sort".
#' @param block Integer or NULL. Block size for batch grouping during
#'   dissection (default NULL = automatic).
#' @param cores Integer. Number of cores for parallel processing (default 1).
#' @param ur Logical. Enable unique combination removal for improved feature
#'   recovery (default TRUE).
#' @param verbose Logical (default TRUE).
#'
#' @return SummarizedExperiment with the new corrected assay added.
#'   If HarmonizR returns fewer features than the input, the SE is subsetted
#'   to match and a warning is issued.
#'
#' @examples
#' \dontrun{
#' source("R/Batch_Diagnostics.R")
#' se_corrected <- batch_correct_proteomics(
#'   se         = result$se_proc,
#'   assay_name = "cycloess",
#'   batch_column = "Batch"
#' )
#' SummarizedExperiment::assayNames(se_corrected)
#' # [1] "raw" "log2" "cycloess" "HarmonizR"
#' }
#' @export
batch_correct_proteomics <- function(
    se,
    assay_name,
    batch_column           = "Batch",
    corrected_assay_name   = "HarmonizR",
    algorithm              = "ComBat",
    ComBat_mode            = 1,
    sort_method            = "sparsity_sort",
    block                  = NULL,
    cores                  = 1,
    ur                     = TRUE,
    verbose                = TRUE
) {

  # --- Check dependencies ---
  .bc_check_harmonizr()

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
    cat("\n=== BATCH CORRECTION (HarmonizR) ===\n")
    cat("- Input assay:", assay_name, "\n")
    cat("- Batch column:", batch_column,
        "(", n_batches, "batches )\n")
    cat("- Algorithm:", algorithm)
    if (algorithm == "ComBat") cat(" (mode", ComBat_mode, ")")
    cat("\n")
  }

  # --- Extract matrix ---
  mat <- SummarizedExperiment::assay(se, assay_name)
  n_features_in <- nrow(mat)
  if (verbose) cat("- Input features:", n_features_in,
                   "| Samples:", ncol(mat), "\n")

  # --- Build description ---
  description <- .bc_build_description(se, batch_column)

  # --- Run HarmonizR ---
  if (verbose) cat("- Running HarmonizR (sort:", sort_method, ") ...\n")
  corrected_mat <- .bc_run_harmonizr(
    mat         = mat,
    description = description,
    algorithm   = algorithm,
    ComBat_mode = ComBat_mode,
    sort_method = sort_method,
    block       = block,
    cores       = cores,
    ur          = ur
  )

  # --- Align output to SE ---
  n_features_out <- nrow(corrected_mat)

  # Handle potential feature loss
  if (n_features_out < n_features_in) {
    n_dropped <- n_features_in - n_features_out
    warning("HarmonizR returned ", n_features_out, " features (", n_dropped,
            " dropped due to batch-specific missingness).",
            "\n  The SE will be subsetted to match.")

    # Subset SE to keep only features returned by HarmonizR
    keep_features <- rownames(corrected_mat)
    se <- se[keep_features, ]
    if (verbose) cat("- Features after correction:", n_features_out,
                     " (", n_dropped, " dropped)\n")
  } else {
    if (verbose) cat("- Features after correction:", n_features_out,
                     " (none dropped)\n")
  }

  # --- Add corrected assay to SE ---
  SummarizedExperiment::assay(se, corrected_assay_name) <- corrected_mat

  if (verbose) {
    cat("- New assay added: '", corrected_assay_name, "'\n", sep = "")
    cat("- Assays in SE:", paste(SummarizedExperiment::assayNames(se),
                                 collapse = ", "), "\n")
    cat("=== BATCH CORRECTION COMPLETADA ===\n")
  }

  se
}
