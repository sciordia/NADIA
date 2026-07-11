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
                           block = NULL,
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
    # covariate can be a single factor or a data.frame of factors
    if (is.data.frame(covariate)) {
      df <- data.frame(condition = condition, covariate)
    } else {
      df <- data.frame(condition = condition, covariate = covariate)
    }
    design <- model.matrix(~ 0 + ., data = df)
    colnames(design)[seq_along(levels(condition))] <- levels(condition)
  }

  # Create contrast matrix
  contrast_strings <- as.character(comparisons)
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_strings,
    levels = design
  )

  # Blocking via duplicateCorrelation
  consensus_cor <- NULL
  if (!is.null(block)) {
    corfit <- tryCatch(
      limma::duplicateCorrelation(data, design, block = block),
      error = function(e) {
        warning("duplicateCorrelation failed: ", conditionMessage(e),
                ". Proceeding without blocking.", call. = FALSE)
        NULL
      }
    )
    if (!is.null(corfit)) consensus_cor <- corfit$consensus.correlation
  }

  # Fit model
  if (!is.null(consensus_cor)) {
    fit <- limma::lmFit(data, design, block = block,
                         correlation = consensus_cor)
  } else {
    fit <- limma::lmFit(data, design)
  }
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
                               block = NULL,
                               eBayes_trend = FALSE, eBayes_robust = FALSE) {
  if (!requireNamespace("limpa", quietly = TRUE)) {
    stop("Para de_method='limpa' necesitas 'limpa'.\n",
         "  BiocManager::install('limpa')")
  }

  if (!is.null(block)) {
    warning("limpa (dpcDE) does not support blocking via duplicateCorrelation. ",
            "bio_replicate_column will be ignored for de_method='limpa'.", call. = FALSE)
  }

  condition <- factor(condition_vector)

  # Design matrix (mismo patron que limma)
  if (is.null(covariate)) {
    design <- model.matrix(~ 0 + condition)
    colnames(design) <- levels(condition)
  } else {
    if (is.data.frame(covariate)) {
      df <- data.frame(condition = condition, covariate)
    } else {
      df <- data.frame(condition = condition, covariate = covariate)
    }
    design <- model.matrix(~ 0 + ., data = df)
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
    # Guarda por signo: con logFC_up = logFC_down = 0, una proteina con
    # logFC == 0 cumpliria >= 0 y <= 0 (marcada Up y luego sobrescrita a Down).
    # Exigir signo estricto la deja como "No Change" (sin direccion).
    df$Change[df$logFC > 0 & df$logFC >= logFC_up   & df[[p_col]] < alpha] <- "Up"
    df$Change[df$logFC < 0 & df$logFC <= logFC_down & df[[p_col]] < alpha] <- "Down"
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
#' @param de_method DE method: "limma" or "limpa" (default: "limma")
#' @param covariate_column Column name(s) in colData for paired/blocked design.
#'   Single string or character vector for multiple covariates (default: NULL)
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
    covariate_column = NULL,
    bio_replicate_column = NULL
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

  # makeContrasts evalua las comparaciones como expresiones (p.ej. "Trt-Ctrl");
  # si un nombre de condicion contiene '-', espacios u otros caracteres no
  # sintacticos (o empieza por digito) el contraste se interpreta mal o falla.
  # Validar temprano con un mensaje claro.
  cond_levels <- unique(as.character(condition_vec))
  bad_levels  <- cond_levels[cond_levels != make.names(cond_levels)]
  if (length(bad_levels) > 0) {
    stop("Nombres de condicion no validos para makeContrasts (contienen '-', ",
         "espacios, u otros caracteres no sintacticos, o empiezan por digito): ",
         paste(bad_levels, collapse = ", "),
         ". Renombralos (p.ej. con make.names) antes del analisis diferencial.")
  }

  # Covariate extraction (supports single or multiple columns)
  covariate <- NULL
  if (!is.null(covariate_column)) {
    missing_cols <- setdiff(covariate_column, names(cd))
    if (length(missing_cols) > 0) {
      stop("Columna(s) de covariable no encontrada(s) en colData del SE: ",
           paste(missing_cols, collapse = ", "))
    }
    # model.matrix hace na.omit por defecto: un NA en la covariable dejaria el
    # design con menos filas que columnas tiene la matriz -> lmFit aborta con
    # un error de dimension poco informativo. Validar explicitamente.
    if (anyNA(cd[, covariate_column, drop = FALSE])) {
      stop("La(s) covariable(s) '", paste(covariate_column, collapse = ", "),
           "' contienen NA en colData. Elimina o imputa esos valores antes del ",
           "analisis diferencial (model.matrix las descartaria y lmFit fallaria).")
    }
    if (length(covariate_column) == 1) {
      covariate <- factor(cd[[covariate_column]])
    } else {
      covariate <- as.data.frame(lapply(cd[covariate_column], factor))
    }
  }

  # Block extraction for duplicateCorrelation
  block <- NULL
  if (!is.null(bio_replicate_column)) {
    if (!bio_replicate_column %in% names(cd))
      stop("bio_replicate_column '", bio_replicate_column, "' not found in colData")
    if (anyNA(cd[[bio_replicate_column]]))
      stop("bio_replicate_column '", bio_replicate_column, "' contiene NA en ",
           "colData; elimina o imputa esos valores antes del analisis diferencial.")
    block_vec <- cd[[bio_replicate_column]]
    if (length(unique(block_vec)) < length(block_vec)) {
      block <- factor(block_vec)
      message("duplicateCorrelation: blocking by '", bio_replicate_column,
              "' (", length(unique(block)), " unique blocks, ", length(block), " samples)")
    } else {
      message("bio_replicate_column '", bio_replicate_column,
              "' \u2014 all values unique, skipping blocking")
    }
  }

  # Run DE analysis
  if (de_method == "limpa") {
    elist <- S4Vectors::metadata(se)$limpa_elist
    if (is.null(elist)) {
      stop("de_method='limpa' requiere imp_method='limpa'. ",
           "No se encontro limpa_elist en metadata del SE.")
    }
    fit <- .perform_limpa_de(elist, condition_vec, comparisons, covariate = covariate,
                              block = block,
                              eBayes_trend = eBayes_trend, eBayes_robust = eBayes_robust)
  } else {
    fit <- .perform_limma(x, condition_vec, comparisons, covariate = covariate,
                          block = block,
                          eBayes_trend = eBayes_trend, eBayes_robust = eBayes_robust)
  }

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
    names(gene_map) <- c("Protein.IDs", "Gene.Names")
    # Deduplicar por Protein.IDs (un mismo ID con dos Gene.Names distintos
    # duplicaria filas del resultado de DE en el merge). Se conserva el primero.
    gene_map <- gene_map[!duplicated(gene_map$Protein.IDs), , drop = FALSE]
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
#' @param de_method DE method: "limma" (default) or "limpa" (probabilistic, requires imp_method="limpa")
#' @param covariate_column Column name(s) in colData for paired/blocked design.
#'   Single string (e.g., "Subject") or character vector (e.g., c("Subject", "Batch")). Default: NULL
#' @param bio_replicate_column Column name in colData identifying biological replicates
#'   (e.g., "Patient", "Subject"). Used with limma::duplicateCorrelation() to account for
#'   technical replicates or paired designs via random effect blocking. Default: NULL
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
    eBayes_trend = NULL,
    eBayes_robust = NULL,
    de_method = "limma",
    covariate_column = NULL,
    bio_replicate_column = NULL,
    condition_column = "Condition",
    verbose = TRUE
) {
  # Validate SE
  stopifnot(inherits(se, "SummarizedExperiment"))

  # Validate de_method
  de_method <- match.arg(de_method, c("limma", "limpa"))

  # Defaults de eBayes segun de_method (NULL = sin fijar por el usuario):
  # limpa usa trend/robust = FALSE (vooma ya modela la tendencia y
  # voomaLmFitWithImputation maneja las proteinas imputadas); limma usa TRUE.
  # Se resuelve aqui (no con missing()) para que funcione tambien cuando el
  # pipeline pasa los argumentos explicitamente.
  default_eb <- de_method != "limpa"
  if (is.null(eBayes_trend))  eBayes_trend  <- default_eb
  if (is.null(eBayes_robust)) eBayes_robust <- default_eb

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
    covariate_column = covariate_column,
    bio_replicate_column = bio_replicate_column
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
