# =============================================================================
# PVCA — Principal Variance Component Analysis
# =============================================================================
#
# Decomposes total variance across experimental factors to identify the main
# sources of variation (batch effects, biological covariates, etc.).
#
# Public functions:
#   pvca_compute()  : Core variance decomposition (returns data.frame)
#   pvca_plot()     : Bar plot of variance components (returns ggplot2)
#   pvca_analysis() : Orchestrator (compute + plot + export)
#
# Based on the PVCA approach described in:
#   - Li et al. (2009) Biostatistics 10(2):317-326
#   - Čuklina et al. (2021) Molecular & Cellular Proteomics 20 (proBatch)
#
# Dependencies (required):
#   - SummarizedExperiment, ggplot2
#
# Dependencies (optional):
#   - pvca   (Bioconductor): pvcaBatchAssess()
#   - Biobase (Bioconductor): ExpressionSet construction
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
#' @param se SummarizedExperiment
#' @param factors Character vector of column names
#' @return invisible(TRUE) or stops with error
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
#'
#' @param se SummarizedExperiment
#' @param assay_name Character, assay to extract
#' @param na_action One of "complete" (drop rows with any NA),
#'   "fill" (replace NAs with fill_value), or "none" (no action).
#' @param fill_value Numeric, used when na_action = "fill" (default -1).
#' @param verbose Logical
#' @return Numeric matrix (proteins x samples) with no NAs
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

#' Build a Biobase ExpressionSet from matrix and SE colData
#' @param mat Numeric matrix (proteins x samples)
#' @param se SummarizedExperiment (for colData extraction)
#' @param factors Character vector of factor column names
#' @return Biobase::ExpressionSet
#' @keywords internal
.pvca_build_eset <- function(mat, se, factors) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  # Keep only the factors and ensure correct sample order
  annot_df <- cd[colnames(mat), factors, drop = FALSE]
  # Convert all to factor (required by pvcaBatchAssess)
  for (col in names(annot_df)) {
    annot_df[[col]] <- as.factor(annot_df[[col]])
  }
  pheno <- Biobase::AnnotatedDataFrame(data = annot_df)
  Biobase::ExpressionSet(assayData = mat, phenoData = pheno)
}

#' Categorize PVCA labels as technical, biological, interaction, or residual
#' @param pvca_df data.frame with columns label, weights
#' @param technical_factors Character vector
#' @param biological_factors Character vector
#' @param variance_threshold Numeric
#' @return data.frame with added column 'category'
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
#' factor plus their interactions and residual variance.
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
  for (pkg in c("pvca", "Biobase")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for pvca_compute(). ",
           "Install from Bioconductor: BiocManager::install('", pkg, "')")
  }

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

  # --- Build ExpressionSet ---
  eset <- .pvca_build_eset(mat, se, factors)

  # --- Run PVCA ---
  if (verbose) message("  Running pvcaBatchAssess (", length(factors),
                       " factors, pca_threshold=", pca_threshold, ") ...")
  pvca_out <- pvca::pvcaBatchAssess(eset, factors, threshold = pca_threshold)

  # --- Format results ---
  pvca_df <- data.frame(
    label   = pvca_out$label,
    weights = as.vector(pvca_out$dat),
    stringsAsFactors = FALSE
  )

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
    message("  PVCA complete. Top components:")
    top <- pvca_df[order(-pvca_df$weights), ]
    for (i in seq_len(min(5, nrow(top)))) {
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
#' source("R/PVCA.R")
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
              title = paste0("PVCA — ", used_assay, " (", n_proteins,
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
