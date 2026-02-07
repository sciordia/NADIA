# =============================================================================
# Differential Expression Analysis Module
# =============================================================================
#
# Functions for differential expression analysis with limma:
#   - Comparison specification (pairwise or vs control)
#   - limma model fitting
#   - Result extraction with fold-change classification
#
# Dependencies:
#   - SummarizedExperiment, S4Vectors
#   - limma
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# --- Null coalescing operator ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# =============================================================================
# INTERNAL FUNCTIONS
# =============================================================================

#' Generate comparisons for differential analysis
#'
#' @param se SummarizedExperiment
#' @param condition_column Condition column (if NULL, uses SE metadata)
#' @param control Control condition. If NULL, generates all pairwise comparisons
#' @return Factor with comparisons in "Treatment-Control" format
#' @keywords internal
.specify_comparisons <- function(
    se,
    condition_column = NULL,
    control = NULL
) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  # Get condition column
  if (is.null(condition_column)) {
    condition_column <- S4Vectors::metadata(se)$condition %||% "Condition"
  }

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if (!condition_column %in% names(cd)) {
    stop("Columna '", condition_column, "' no encontrada en colData")
  }

  conditions <- unique(as.character(cd[[condition_column]]))

  if (!is.null(control)) {
    if (!control %in% conditions) {
      stop("Control '", control, "' no esta en las condiciones: ",
           paste(conditions, collapse = ", "))
    }
    # Comparisons vs control
    others <- setdiff(conditions, control)
    comparisons <- paste0(others, "-", control)
  } else {
    # All pairwise comparisons
    comparisons <- character()
    for (i in seq_along(conditions)) {
      for (j in seq_along(conditions)) {
        if (i < j) {
          comparisons <- c(comparisons, paste0(conditions[j], "-", conditions[i]))
        }
      }
    }
  }

  factor(comparisons)
}

#' Run limma analysis
#'
#' @param data Log2 intensity matrix
#' @param condition_vector Condition vector aligned with columns
#' @param comparisons Comparison vector
#' @param covariate Optional covariate for the model
#' @param eBayes_trend Use trend estimation in eBayes (default: TRUE)
#' @param eBayes_robust Use robust estimation in eBayes (default: TRUE)
#' @return limma fit object
#' @keywords internal
.perform_limma <- function(data, condition_vector, comparisons, covariate = NULL,
                           eBayes_trend = TRUE, eBayes_robust = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Se requiere el paquete 'limma'")
  }

  condition <- factor(condition_vector)

  # Create design matrix
  if (is.null(covariate)) {
    design <- model.matrix(~ 0 + condition)
    colnames(design) <- levels(condition)
  } else {
    design <- model.matrix(~ 0 + condition + covariate)
    colnames(design)[seq_along(levels(condition))] <- levels(condition)
  }

  # Create contrast matrix
  contrast_strings <- as.character(comparisons)
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_strings,
    levels = design
  )

  # Fit model
  fit <- limma::lmFit(data, design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit, trend = eBayes_trend, robust = eBayes_robust)

  fit
}

#' Extract results from limma fit
#'
#' @param fit limma fit object
#' @param comparisons Comparison vector
#' @param logFC_up Upper logFC threshold for "Up"
#' @param logFC_down Lower logFC threshold for "Down"
#' @param alpha Significance threshold
#' @param p_adj Use adjusted p-value (TRUE) or raw (FALSE)
#' @return Data frame with results
#' @keywords internal
.extract_limma_results <- function(
    fit,
    comparisons,
    logFC_up = 1,
    logFC_down = -1,
    alpha = 0.05,
    p_adj = TRUE
) {
  results_list <- lapply(seq_along(comparisons), function(i) {
    comp <- as.character(comparisons[i])
    tt <- limma::topTable(fit, coef = i, number = Inf, sort.by = "none")

    # Standard columns
    df <- data.frame(
      Protein.IDs = rownames(tt),
      logFC = tt$logFC,
      P.Value = tt$P.Value,
      adj.P.Val = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )

    # Classify changes
    p_col <- if (p_adj) "adj.P.Val" else "P.Value"
    df$Change <- "No Change"
    df$Change[df$logFC >= logFC_up & df[[p_col]] < alpha] <- "Up"
    df$Change[df$logFC <= logFC_down & df[[p_col]] < alpha] <- "Down"
    df$Change <- factor(df$Change, levels = c("Up", "Down", "No Change"))

    df$Comparison <- comp
    df
  })

  do.call(rbind, results_list)
}

#' Run differential expression analysis
#'
#' Internal function that coordinates DE analysis with limma.
#'
#' @param se SummarizedExperiment with processed data
#' @param comparisons Comparisons to perform (result of .specify_comparisons)
#' @param assay_name Assay name to use. If NULL, uses the last available
#' @param condition_column Condition column. If NULL, uses SE metadata
#' @param logFC Apply logFC filter (default: TRUE)
#' @param logFC_up Upper logFC threshold (default: 1)
#' @param logFC_down Lower logFC threshold (default: -1)
#' @param p_adj Use adjusted p-value (default: TRUE)
#' @param alpha Significance threshold (default: 0.05)
#' @param eBayes_trend Use trend estimation in eBayes (default: TRUE)
#' @param eBayes_robust Use robust estimation in eBayes (default: TRUE)
#' @return Data frame with DE results
#' @keywords internal
.run_DE <- function(
    se,
    comparisons,
    assay_name = NULL,
    condition_column = NULL,
    logFC = TRUE,
    logFC_up = 1,
    logFC_down = -1,
    p_adj = TRUE,
    alpha = 0.05,
    eBayes_trend = TRUE,
    eBayes_robust = TRUE
) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  # Determine assay
  if (is.null(assay_name)) {
    assay_names <- SummarizedExperiment::assayNames(se)
    assay_name <- assay_names[length(assay_names)]
  }

  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", assay_name, "' no encontrado")
  }

  # Get data
  x <- SummarizedExperiment::assay(se, assay_name)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  rd <- as.data.frame(SummarizedExperiment::rowData(se))

  # Condition column
  if (is.null(condition_column)) {
    condition_column <- S4Vectors::metadata(se)$condition %||% "Condition"
  }

  condition_vec <- cd[[condition_column]]

  # Run limma
  fit <- .perform_limma(x, condition_vec, comparisons, covariate = NULL,
                        eBayes_trend = eBayes_trend, eBayes_robust = eBayes_robust)

  # Extract results
  if (!logFC) {
    logFC_up <- 0
    logFC_down <- 0
  }

  results <- .extract_limma_results(
    fit, comparisons,
    logFC_up = logFC_up,
    logFC_down = logFC_down,
    alpha = alpha,
    p_adj = p_adj
  )

  # Add gene information
  if ("Gene.Names" %in% names(rd)) {
    gene_map <- rd[, c("Protein.IDs", "Gene.Names"), drop = FALSE]
    gene_map <- unique(gene_map)
    names(gene_map) <- c("Protein.IDs", "Gene.Names")
    results <- merge(results, gene_map, by = "Protein.IDs", all.x = TRUE, sort = FALSE)
  }

  # Add Assay column
  results$Assay <- assay_name

  # Reorder columns
  col_order <- c("Protein.IDs", "Gene.Names", "logFC", "P.Value",
                 "adj.P.Val", "Change", "Comparison", "Assay")
  col_order <- intersect(col_order, names(results))
  results <- results[, col_order]

  results
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Perform differential expression analysis on proteomics data
#'
#' Coordinates comparison specification, limma analysis, and result extraction.
#'
#' @param se SummarizedExperiment with imputed assay (output of impute_proteomics)
#' @param assay_name Assay name to use. If NULL, uses the last available assay
#' @param comparisons Comparisons to perform. If NULL, auto-generates (pairwise or vs control)
#' @param control Control condition. If NULL, generates all pairwise comparisons
#' @param logFC_threshold LogFC threshold for significance (default: 0)
#' @param alpha Adjusted p-value threshold (default: 0.05)
#' @param p_adj Use adjusted p-value (default: TRUE)
#' @param eBayes_trend Use trend estimation in eBayes (default: TRUE, recommended for proteomics)
#' @param eBayes_robust Use robust estimation in eBayes (default: TRUE, recommended for proteomics)
#' @param condition_column Condition column name (default: "Condition")
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return List with:
#'   \itemize{
#'     \item DEPs_results: Data frame with differential expression results
#'     \item comparisons: Comparisons performed
#'     \item assay_name: Assay name used
#'   }
#'
#' @examples
#' \dontrun{
#' de_result <- de_analysis_proteomics(
#'   se = imp_result$se,
#'   assay_name = "Cycloess",
#'   control = "A",
#'   alpha = 0.05
#' )
#' }
#'
#' @export
de_analysis_proteomics <- function(
    se,
    assay_name = NULL,
    comparisons = NULL,
    control = NULL,
    logFC_threshold = 0,
    alpha = 0.05,
    p_adj = TRUE,
    eBayes_trend = TRUE,
    eBayes_robust = TRUE,
    condition_column = "Condition",
    verbose = TRUE
) {
  # Validate SE
  stopifnot(inherits(se, "SummarizedExperiment"))

  # Validate required packages
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Se requiere el paquete 'limma'. ",
         "Instalalo con BiocManager::install('limma')")
  }

  # Determine assay name
  if (is.null(assay_name)) {
    available_assays <- SummarizedExperiment::assayNames(se)
    assay_name <- available_assays[length(available_assays)]
  }

  if (verbose) cat("\n=== ANALISIS DIFERENCIAL (limma) ===\n")

  # Generate comparisons if not specified
  if (is.null(comparisons)) {
    comparisons <- .specify_comparisons(se, condition_column = condition_column,
                                        control = control)
  }

  if (verbose) cat("- Comparaciones:", paste(comparisons, collapse = ", "), "\n")

  # Run DE analysis
  DEPs_results <- .run_DE(
    se = se,
    comparisons = comparisons,
    assay_name = assay_name,
    condition_column = condition_column,
    logFC = TRUE,
    logFC_up = logFC_threshold,
    logFC_down = -logFC_threshold,
    p_adj = p_adj,
    alpha = alpha,
    eBayes_trend = eBayes_trend,
    eBayes_robust = eBayes_robust
  )

  if (verbose) {
    n_sig <- sum(DEPs_results$Change != "No Change")
    cat("- Proteinas diferenciales (total):", n_sig, "\n")

    for (comp in unique(DEPs_results$Comparison)) {
      subset <- DEPs_results[DEPs_results$Comparison == comp, ]
      n_up <- sum(subset$Change == "Up")
      n_down <- sum(subset$Change == "Down")
      cat("  ", comp, ": Up=", n_up, ", Down=", n_down, "\n", sep = "")
    }
  }

  # Return
  list(
    DEPs_results = DEPs_results,
    comparisons = comparisons,
    assay_name = assay_name
  )
}
