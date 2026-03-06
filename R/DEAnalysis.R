# =============================================================================
# Differential Expression Analysis Module
# =============================================================================
#
# Functions for differential expression analysis with limma/limpa/LimROTS:
#   - Comparison specification (pairwise or vs control)
#   - limma/limpa model fitting, LimROTS bootstrapped reproducibility
#   - Result extraction with fold-change classification
#
# Dependencies:
#   - SummarizedExperiment, S4Vectors
#   - limma
#   - Optional: limpa, LimROTS, BiocParallel
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

#' Run limpa dpcDE analysis
#'
#' Uses limpa::dpcDE with precision weights from the EList (standard errors)
#' produced by dpcQuantByRow. The EList is stored in SE metadata by
#' impute_proteomics() when imp_method="limpa".
#'
#' @param elist EList from limpa (with $E and $weights)
#' @param condition_vector Condition vector aligned with columns
#' @param comparisons Comparison vector
#' @param covariate Optional covariate for the model
#' @param eBayes_trend Use trend estimation in eBayes (default: FALSE, vooma already models the trend)
#' @param eBayes_robust Use robust estimation in eBayes (default: FALSE, voomaLmFitWithImputation
#'   already handles imputed proteins)
#' @return limma MArrayLM fit object
#' @keywords internal
.perform_limpa_de <- function(elist, condition_vector, comparisons, covariate = NULL,
                               eBayes_trend = FALSE, eBayes_robust = FALSE) {
  if (!requireNamespace("limpa", quietly = TRUE)) {
    stop("Para de_method='limpa' necesitas 'limpa'.\n",
         "  BiocManager::install('limpa')")
  }

  condition <- factor(condition_vector)

  # Design matrix (mismo patron que limma)
  if (is.null(covariate)) {
    design <- model.matrix(~ 0 + condition)
    colnames(design) <- levels(condition)
  } else {
    design <- model.matrix(~ 0 + condition + covariate)
    colnames(design)[seq_along(levels(condition))] <- levels(condition)
  }

  # Contrast matrix
  contrast_strings <- as.character(comparisons)
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_strings,
    levels = design
  )

  # dpcDE: ajuste con precision weights de SEs
  # voomaLmFitWithImputation ya modela la tendencia de varianza via vooma,
  # por lo que eBayes defaults son FALSE/FALSE (vignette Li, Cobbold, Smyth 2025).
  # Se exponen como configurables para usuarios avanzados (como hace msdap).
  fit <- limpa::dpcDE(elist, design, plot = FALSE)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit, trend = eBayes_trend, robust = eBayes_robust)

  fit
}

#' Run LimROTS analysis on a 2-group subset
#'
#' LimROTS does not support pairwise contrasts; it requires a 2-group SE.
#' This function runs one comparison at a time on a subsetted SE.
#'
#' @param se_subset SummarizedExperiment with exactly 2 conditions
#' @param condition_column Condition column name in colData
#' @param comparison_label Label for the comparison (e.g. "B-A")
#' @param niter Bootstrap iterations (default: 1000)
#' @param K Top features for reproducibility ranking (default: NULL = nrow/4)
#' @param eBayes_trend Use trend estimation in eBayes (default: TRUE)
#' @param eBayes_robust Use robust estimation in eBayes (default: TRUE)
#' @param BPPARAM BiocParallel param (default: NULL = SerialParam)
#' @param verbose Print progress (default: TRUE)
#' @return SummarizedExperiment with LimROTS results in rowData
#' @keywords internal
.perform_LimROTS_de <- function(se_subset, condition_column, comparison_label,
                                 niter = 1000, K = NULL,
                                 eBayes_trend = TRUE, eBayes_robust = TRUE,
                                 BPPARAM = NULL, verbose = TRUE) {
  if (!requireNamespace("LimROTS", quietly = TRUE)) {
    stop("Para de_method='LimROTS' necesitas 'LimROTS'.\n",
         "  BiocManager::install('LimROTS')")
  }

  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::SerialParam()
  }

  if (is.null(K)) K <- floor(nrow(se_subset) / 4)

  se_result <- LimROTS::LimROTS(
    x               = se_subset,
    niter           = niter,
    K               = K,
    meta.info       = c(condition_column),
    group.name      = condition_column,
    formula.str     = paste0("~ 0 + ", condition_column),
    trend           = eBayes_trend,
    robust          = eBayes_robust,
    BPPARAM         = BPPARAM,
    verbose         = verbose,
    log             = TRUE
  )

  se_result
}

#' Extract results from LimROTS SE output
#'
#' Formats LimROTS rowData into the same structure as .extract_limma_results().
#'
#' @param se_result SummarizedExperiment returned by LimROTS
#' @param comparison_label Label for the comparison
#' @param logFC_up Upper logFC threshold for "Up"
#' @param logFC_down Lower logFC threshold for "Down"
#' @param alpha Significance threshold
#' @param p_adj Use FDR (TRUE) or raw p-value (FALSE)
#' @return Data frame with Protein.IDs, logFC, P.Value, adj.P.Val, Change, Comparison
#' @keywords internal
.extract_LimROTS_results <- function(se_result, comparison_label,
                                      logFC_up = 1, logFC_down = -1,
                                      alpha = 0.05, p_adj = TRUE) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se_result))

  df <- data.frame(
    Protein.IDs = rownames(rd),
    logFC       = rd$corrected.logfc,
    P.Value     = rd$pvalue,
    adj.P.Val   = rd$FDR,
    stringsAsFactors = FALSE
  )

  # Classify changes (same logic as .extract_limma_results)
  p_col <- if (p_adj) "adj.P.Val" else "P.Value"
  df$Change <- "No Change"
  df$Change[df$logFC >= logFC_up & df[[p_col]] < alpha] <- "Up"
  df$Change[df$logFC <= logFC_down & df[[p_col]] < alpha] <- "Down"
  df$Change <- factor(df$Change, levels = c("Up", "Down", "No Change"))

  df$Comparison <- comparison_label
  df
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
#' @param de_method DE method: "limma", "limpa", or "LimROTS" (default: "limma")
#' @param niter Bootstrap iterations for LimROTS (default: 1000)
#' @param K Top features for LimROTS reproducibility ranking (default: NULL = nrow/4)
#' @param BPPARAM BiocParallel param for LimROTS (default: NULL = serial)
#' @param verbose Print progress messages (default: TRUE)
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
    eBayes_robust = TRUE,
    de_method = "limma",
    niter = 1000,
    K = NULL,
    BPPARAM = NULL,
    verbose = TRUE
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

  # Apply logFC thresholds
  if (!logFC) {
    logFC_up <- 0
    logFC_down <- 0
  }

  # Run DE analysis
  if (de_method == "LimROTS") {
    # LimROTS: run separately for each 2-group comparison
    all_conditions <- unique(as.character(condition_vec))

    if (verbose) {
      cat("  LimROTS:", length(comparisons), "comparaciones x", niter, "bootstraps\n")
    }

    results_list <- lapply(seq_along(comparisons), function(i) {
      comp <- as.character(comparisons[i])

      # Parse comparison "B-A" -> treatment="B", control="A"
      cond_treatment <- NULL
      cond_control <- NULL
      for (cond in all_conditions) {
        if (startsWith(comp, paste0(cond, "-"))) {
          cond_treatment <- cond
          cond_control <- sub(paste0("^", cond, "-"), "", comp)
          break
        }
      }
      if (is.null(cond_treatment)) stop("No se pudo parsear la comparacion: ", comp)

      # Subset SE to the 2 conditions
      keep_samples <- condition_vec %in% c(cond_treatment, cond_control)
      se_2group <- se[, keep_samples]

      # Filter rows with any NA in the assay (avoid NA coefficients in bootstrap)
      assay_mat <- SummarizedExperiment::assay(se_2group, assay_name)
      complete_rows <- rowSums(is.na(assay_mat)) == 0
      if (sum(!complete_rows) > 0 && verbose) {
        cat("    (filtrando", sum(!complete_rows), "proteinas con NAs)\n")
      }
      se_2group <- se_2group[complete_rows, ]

      # Set factor levels: treatment FIRST (LimROTS: group1 - group2)
      cd_2group <- as.data.frame(SummarizedExperiment::colData(se_2group))
      cd_2group[[condition_column]] <- factor(
        cd_2group[[condition_column]],
        levels = c(cond_treatment, cond_control)
      )
      SummarizedExperiment::colData(se_2group)[[condition_column]] <- cd_2group[[condition_column]]

      if (verbose) cat("  - Ejecutando LimROTS para", comp, "\n")

      # Run LimROTS
      se_result <- .perform_LimROTS_de(
        se_subset        = se_2group,
        condition_column = condition_column,
        comparison_label = comp,
        niter            = niter,
        K                = K,
        eBayes_trend     = eBayes_trend,
        eBayes_robust    = eBayes_robust,
        BPPARAM          = BPPARAM,
        verbose          = verbose
      )

      # Extract results in standard format
      .extract_LimROTS_results(se_result, comp,
        logFC_up = logFC_up, logFC_down = logFC_down,
        alpha = alpha, p_adj = p_adj)
    })

    results <- do.call(rbind, results_list)

  } else {
    # limma or limpa
    if (de_method == "limpa") {
      elist <- S4Vectors::metadata(se)$limpa_elist
      if (is.null(elist)) {
        stop("de_method='limpa' requiere imp_method='limpa'. ",
             "No se encontro limpa_elist en metadata del SE.")
      }
      fit <- .perform_limpa_de(elist, condition_vec, comparisons, covariate = NULL,
                                eBayes_trend = eBayes_trend, eBayes_robust = eBayes_robust)
    } else {
      fit <- .perform_limma(x, condition_vec, comparisons, covariate = NULL,
                            eBayes_trend = eBayes_trend, eBayes_robust = eBayes_robust)
    }

    results <- .extract_limma_results(
      fit, comparisons,
      logFC_up = logFC_up,
      logFC_down = logFC_down,
      alpha = alpha,
      p_adj = p_adj
    )
  }

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
#' @param de_method DE method: "limma" (default), "limpa" (probabilistic, requires
#'   imp_method="limpa"), or "LimROTS" (bootstrapped reproducibility-optimized test statistic)
#' @param LimROTS_niter Bootstrap iterations for LimROTS (default: 1000). Higher = more
#'   precise but slower.
#' @param LimROTS_K Top features for reproducibility ranking (default: NULL = nrow/4)
#' @param LimROTS_BPPARAM BiocParallel param for LimROTS (default: NULL = serial).
#'   Use BiocParallel::MulticoreParam(4) for parallel.
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
#'   assay_name = "ImpSeqRob_Min",
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
    de_method = "limma",
    LimROTS_niter   = 1000,
    LimROTS_K       = NULL,
    LimROTS_BPPARAM = NULL,
    condition_column = "Condition",
    verbose = TRUE
) {
  # Validate SE
  stopifnot(inherits(se, "SummarizedExperiment"))

  # Validate de_method
  de_method <- match.arg(de_method, c("limma", "limpa", "LimROTS"))

  # Para limpa, defaults de eBayes son FALSE (vooma ya modela la tendencia)
  if (de_method == "limpa") {
    if (missing(eBayes_trend))  eBayes_trend  <- FALSE
    if (missing(eBayes_robust)) eBayes_robust <- FALSE
  }

  # Para LimROTS, trend=FALSE por defecto: el bootstrap interno puede generar
  # fits parciales con NA coefficients, y trend=TRUE usa Amean como covariable
  # que hereda esos NAs, crasheando fitFDistUnequalDF1.
  if (de_method == "LimROTS") {
    if (missing(eBayes_trend))  eBayes_trend  <- FALSE
    if (missing(eBayes_robust)) eBayes_robust <- TRUE
  }

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

  if (verbose) cat("\n=== ANALISIS DIFERENCIAL (", de_method, ") ===\n")

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
    eBayes_robust = eBayes_robust,
    de_method = de_method,
    niter   = LimROTS_niter,
    K       = LimROTS_K,
    BPPARAM = LimROTS_BPPARAM,
    verbose = verbose
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
