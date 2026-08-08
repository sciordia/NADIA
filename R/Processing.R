# =============================================================================
# Proteomics Data Processing (Coordinator)
# =============================================================================
#
# Main function that coordinates the complete pipeline:
#   1. Normalization.R - Filtering, multiple normalization methods
#   2. Imputation.R   - 17 imputation methods (combo MAR+MNAR + individual)
#   3. DEAnalysis.R   - Differential analysis with limma
#
# Dependencies: see the individual modules
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================

# Normalization.R, Imputation.R, DEAnalysis.R and Batch_Correction.R share a
# namespace with this file, so there is no need to load them.

# =============================================================================
# LINKER FUNCTIONS (internal)
# =============================================================================

#' Prepare metadata from proteomics_data object
#'
#' @param preprocessing proteomics_data list (output of preprocess_spectronaut or preprocess_tmt)
#' @return Data frame with columns: Column, Condition, Replicate
#' @keywords internal
.prepare_metadata <- function(preprocessing, covariate_df = NULL) {
  stopifnot(inherits(preprocessing, "proteomics_data"))

  md <- preprocessing$metadata

  result <- data.frame(
    Column = md$Coding,
    Condition = md$R.Condition,
    Replicate = md$R.Replicate,
    stringsAsFactors = FALSE
  )

  # Merge covariate information if provided

  if (!is.null(covariate_df)) {
    if (!"Column" %in% names(covariate_df)) {
      stop("covariate_df must contain a 'Column' column")
    }
    # Ignore columns already derived from the preprocessing (Condition, Replicate, ...)
    redundant <- setdiff(intersect(names(covariate_df), names(result)), "Column")
    if (length(redundant) > 0) {
      message("Ignoring covariate_df columns already present in metadata: ",
              paste(redundant, collapse = ", "))
      covariate_df <- covariate_df[, setdiff(names(covariate_df), redundant),
                                   drop = FALSE]
    }
    if (ncol(covariate_df) > 1) {
      result <- merge(result, covariate_df, by = "Column", all.x = TRUE)
      # Check for unmatched samples
      new_cols <- setdiff(names(covariate_df), "Column")
      na_check <- vapply(new_cols, function(col) any(is.na(result[[col]])),
                         logical(1))
      if (any(na_check)) {
        missing_cols <- names(na_check)[na_check]
        stop("covariate_df does not cover all samples. NAs in: ",
             paste(missing_cols, collapse = ", "))
      }
    }
  }

  rownames(result) <- result$Column
  result
}

#' Prepare protein data from proteomics_data object
#'
#' @param preprocessing proteomics_data list (output of preprocess_spectronaut or preprocess_tmt)
#' @return Data frame with columns: ProteinGroups, GeneNames, UniqPepts, and intensity columns
#' @keywords internal
.prepare_protein_data <- function(preprocessing) {
  stopifnot(inherits(preprocessing, "proteomics_data"))

  pq <- preprocessing$protein_quant

  # Identify quantity columns (PG.Quantity_*)
  quantity_cols <- grep("^PG\\.Quantity_", names(pq), value = TRUE)
  if (length(quantity_cols) == 0) {
    stop("No PG.Quantity_* columns found in protein_quant")
  }

  # Identify unique peptide columns
  pept_cols <- grep("^PG\\.NrOfStrippedSequencesUsedForQuantification_", names(pq), value = TRUE)

  # Compute UniqPepts as max per row
  if (length(pept_cols) > 0) {
    pept_mat <- as.matrix(pq[, pept_cols, drop = FALSE])
    UniqPepts <- apply(pept_mat, 1, function(x) max(x, na.rm = TRUE))
    UniqPepts[!is.finite(UniqPepts)] <- 0L
  } else {
    UniqPepts <- rep(NA_integer_, nrow(pq))
  }

  # Build result
  result <- data.frame(
    ProteinGroups = pq$PG.ProteinGroups,
    GeneNames = pq$PG.Genes,
    UniqPepts = as.integer(UniqPepts),
    stringsAsFactors = FALSE
  )

  # Add intensity columns with simplified names (Condition_Replicate)
  intensity_data <- pq[, quantity_cols, drop = FALSE]
  names(intensity_data) <- gsub("^PG\\.Quantity_", "", names(intensity_data))

  result <- cbind(result, intensity_data)
  result
}

# =============================================================================
# EXPORT UTILITY FUNCTIONS (internal)
# =============================================================================

#' Convert SummarizedExperiment to long format
#'
#' @param se SummarizedExperiment
#' @param assay_names Assay names to include. If NULL, uses all
#' @return Data frame in long format
#' @keywords internal
.se_to_long <- function(se, assay_names = NULL) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  available_assays <- SummarizedExperiment::assayNames(se)
  if (is.null(assay_names)) {
    assay_names <- available_assays
  } else {
    assay_names <- intersect(assay_names, available_assays)
    if (length(assay_names) == 0) {
      stop("None of the specified assays is available")
    }
  }

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  rd <- as.data.frame(SummarizedExperiment::rowData(se))

  result_list <- lapply(assay_names, function(assay_name) {
    mat <- SummarizedExperiment::assay(se, assay_name)

    df_list <- lapply(seq_len(ncol(mat)), function(j) {
      sample_name <- colnames(mat)[j]
      data.frame(
        Column = sample_name,
        Assay = assay_name,
        Intensity = mat[, j],
        Protein.IDs = rownames(mat),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, df_list)
  })

  long_df <- do.call(rbind, result_list)

  if ("Condition" %in% names(cd)) {
    long_df$Condition <- cd[long_df$Column, "Condition"]
  }
  if ("Replicate" %in% names(cd)) {
    long_df$Replicate <- cd[long_df$Column, "Replicate"]
  }

  col_order <- c("Column", "Assay", "Intensity", "Condition", "Replicate", "Protein.IDs")
  col_order <- intersect(col_order, names(long_df))
  long_df <- long_df[, col_order, drop = FALSE]

  rownames(long_df) <- NULL
  long_df
}

#' Prepare PCA input with differential expression information
#'
#' @param se SummarizedExperiment
#' @param DEPs_results Data frame with DE results
#' @param assay_name Assay name to use
#' @param alpha Significance threshold for sig_any
#' @return Data frame with intensities and adjP columns per comparison
#' @keywords internal
.prepare_pca_input <- function(se, DEPs_results, assay_name, alpha = 0.05) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", assay_name, "' not found in SE")
  }

  mat <- SummarizedExperiment::assay(se, assay_name)
  cd <- as.data.frame(SummarizedExperiment::colData(se))

  long_list <- lapply(seq_len(ncol(mat)), function(j) {
    sample_name <- colnames(mat)[j]
    data.frame(
      SampleID = sample_name,
      FeatureID = rownames(mat),
      Intensity = mat[, j],
      stringsAsFactors = FALSE
    )
  })
  long_df <- do.call(rbind, long_list)

  if ("Condition" %in% names(cd)) {
    long_df$Condition <- cd[long_df$SampleID, "Condition"]
  }
  if ("Replicate" %in% names(cd)) {
    long_df$Replicate <- cd[long_df$SampleID, "Replicate"]
  }

  comparisons <- unique(DEPs_results$Comparison)

  for (comp in comparisons) {
    subset_de <- DEPs_results[DEPs_results$Comparison == comp, ]
    adjp_map <- setNames(subset_de$adj.P.Val, subset_de$Protein.IDs)
    col_name <- paste0("adjP_", comp)
    long_df[[col_name]] <- adjp_map[long_df$FeatureID]
  }

  adjp_cols <- grep("^adjP_", names(long_df), value = TRUE)
  if (length(adjp_cols) > 0) {
    adjp_mat <- as.matrix(long_df[, adjp_cols, drop = FALSE])
    long_df$sig_any <- apply(adjp_mat, 1, function(x) any(x < alpha, na.rm = TRUE))
  } else {
    long_df$sig_any <- FALSE
  }

  rownames(long_df) <- NULL
  long_df
}

#' Export data frame to TSV and/or Parquet
#'
#' @param data Data frame to export
#' @param filepath_base Base path without extension
#' @param format "tsv", "parquet", or "both"
#' @return Invisible NULL
#' @keywords internal
.export_data <- function(data, filepath_base, format = "tsv") {
  format <- match.arg(format, c("tsv", "parquet", "both"))

  if (format %in% c("tsv", "both")) {
    tsv_file <- paste0(filepath_base, ".tsv")
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_tsv(data, tsv_file)
    } else {
      write.table(data, tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  if (format %in% c("parquet", "both")) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      warning("Package 'arrow' is not installed. Cannot export to Parquet. ",
              "Install it with: install.packages('arrow')")
    } else {
      parquet_file <- paste0(filepath_base, ".parquet")
      arrow::write_parquet(data, parquet_file)
    }
  }

  invisible(NULL)
}

#' Compute missing-value percentages per protein
#'
#' Four percentages, all measured on the normalized assay, that is, before
#' imputation, and all going from the widest scope to the narrowest:
#' - MissGlobal: % NAs over every sample in the experiment, including the
#'   conditions that take no part in the comparison. It does not depend on the
#'   comparison, so a given protein carries the same value in all of them.
#' - MissComp: % NAs over the replicates of the two conditions being compared
#'   only. With more than two conditions this is not the same as MissGlobal.
#' - MissCND1: % NAs over the samples of Cond1 (numerator, before "-")
#' - MissCND2: % NAs over the samples of Cond2 (denominator, after "-")
#'
#' @param se SummarizedExperiment with normalized assay
#' @param DEPs_results Data frame with DE results (must have Protein.IDs, Comparison)
#' @param norm_assay_name Name of the normalized assay
#' @return DEPs_results with MissGlobal, MissComp, MissCND1, MissCND2 columns added
#' @keywords internal
.compute_missing_pct <- function(se, DEPs_results, norm_assay_name) {
  x_norm <- SummarizedExperiment::assay(se, norm_assay_name)

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cond_samples <- split(cd$Column, cd$Condition)

  # Global missingness: one value per protein over every sample in the
  # experiment. It does not depend on the comparison, so it is computed once
  # here and joined by protein alone, rather than repeated inside the per
  # comparison loop below.
  na_global <- if (ncol(x_norm) > 0) {
    round(100 * rowSums(is.na(x_norm)) / ncol(x_norm), 2)
  } else {
    rep(NA_real_, nrow(x_norm))
  }

  known <- names(cond_samples)
  comps <- unique(as.character(DEPs_results$Comparison))
  pct_list <- lapply(comps, function(comp) {
    # comp = "cond1-cond2" (numerator-denominator). Derive cond1/cond2 by matching
    # against the known condition names instead of blindly re-parsing on "-"
    # (robust even if a name were to contain a hyphen).
    cond1 <- NA_character_
    cond2 <- NA_character_
    for (c1 in known) {
      prefix <- paste0(c1, "-")
      if (startsWith(comp, prefix)) {
        rest <- substring(comp, nchar(prefix) + 1L)
        if (rest %in% known) { cond1 <- c1; cond2 <- rest; break }
      }
    }
    if (is.na(cond1)) {  # fallback: plain split on "-"
      parts <- trimws(strsplit(comp, "-")[[1]])
      cond1 <- parts[1]  # numerator (e.g. B in "B-A")
      cond2 <- parts[2]  # denominator (e.g. A in "B-A")
    }

    s1 <- intersect(cond_samples[[cond1]] %||% character(0), colnames(x_norm))
    s2 <- intersect(cond_samples[[cond2]] %||% character(0), colnames(x_norm))
    s_all <- c(s1, s2)

    n1 <- length(s1); n2 <- length(s2); n_all <- length(s_all)

    if (n_all == 0) {
      return(data.frame(Protein.IDs = rownames(x_norm), Comparison = comp,
                        MissComp = NA_real_, MissCND1 = NA_real_,
                        MissCND2 = NA_real_, stringsAsFactors = FALSE))
    }

    na_all <- if (n_all > 0) rowSums(is.na(x_norm[, s_all, drop = FALSE])) else rep(0L, nrow(x_norm))
    na1    <- if (n1 > 0)    rowSums(is.na(x_norm[, s1, drop = FALSE]))    else rep(NA_real_, nrow(x_norm))
    na2    <- if (n2 > 0)    rowSums(is.na(x_norm[, s2, drop = FALSE]))    else rep(NA_real_, nrow(x_norm))

    data.frame(
      Protein.IDs = rownames(x_norm),
      Comparison  = comp,
      MissComp    = round(100 * na_all / n_all, 2),
      MissCND1    = if (n1 > 0) round(100 * na1 / n1, 2) else NA_real_,
      MissCND2    = if (n2 > 0) round(100 * na2 / n2, 2) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  pct_df <- do.call(rbind, pct_list)

  DEPs_results <- merge(DEPs_results, pct_df,
                        by = c("Protein.IDs", "Comparison"),
                        all.x = TRUE, sort = FALSE)

  # match(), not %in%: Protein.IDs repeats across comparisons in DEPs_results
  DEPs_results$MissGlobal <-
    na_global[match(DEPs_results$Protein.IDs, rownames(x_norm))]

  missing_cols <- c("MissGlobal", "MissComp", "MissCND1", "MissCND2")
  other_cols <- setdiff(names(DEPs_results), missing_cols)
  DEPs_results[, c(other_cols, intersect(missing_cols, names(DEPs_results)))]
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Process proteomics data from Spectronaut
#'
#' Complete processing pipeline: filtering, normalization, imputation
#' and differential expression analysis. Coordinates sub-modules
#' Normalization.R, Imputation.R, and DEAnalysis.R.
#'
#' @param preprocessing proteomics_data list (result of preprocess_spectronaut or preprocess_tmt)
#' @param export_dir Output directory for exported files. Defaults to `NULL`,
#'   which writes nothing to disk; pass a path to enable the exports controlled
#'   by `export_normalized`, `export_imputed`, `export_volcano`,
#'   `export_boxplot` and `export_pca`.
#' @param min_reps_filter Minimum replicates for filtering. If NULL, auto-computed
#' @param min_groups_filter Minimum groups for filtering (default: 1)
#' @param norm_method Normalization method passed to normalize_proteomics()
#'   (default: "cycloess"). See normalize_proteomics() for all 14 options.
#' @param cyclic_loess_method Cyclic Loess method: "fast" or "pairs" (default: "fast")
#' @param cyclic_loess_iterations Number of iterations for Cyclic Loess (default: 3)
#' @param cyclic_loess_span Span parameter for Cyclic Loess (default: 0.7)
#' @param batch_correct Logical. Apply BERT batch correction after
#'   normalization (default: FALSE). Requires covariate_df with batch_column.
#' @param batch_column Column in covariate_df containing batch assignments
#'   (default: "Batch"). Must have >= 2 unique values.
#' @param batch_algorithm Batch correction algorithm: "ComBat" (default), "limma", or "ref"
#' @param batch_ComBat_mode Integer 1-4 for ComBat parametric/mean-only settings (default: 1)
#' @param batch_covariates Character vector of colData column names to use as
#'   categorical covariates for BERT batch correction (default: NULL)
#' @param batch_qualitycontrol Logical. Compute ASW quality metrics (default: FALSE)
#' @param imp_method Imputation method (default: "combo"). See impute_proteomics() for all options.
#' @param mar_method MAR method for combo mode (default: "Impseqrob")
#' @param mnar_method MNAR method for combo mode (default: "min")
#' @param prop_na_mnar NA proportion for MNAR classification (default: 0.51)
#' @param prop_present_mar Present proportion for MAR (default: 0.5)
#' @param min_present_mar Minimum present values for MAR (default: 1)
#' @param require_n_conditions Required conditions with presence (default: 1)
#' @param max_na_prop Maximum NA proportion above which a protein is removed
#'   before imputation. `NULL` (the default) disables the filter.
#' @param method_args Named list of per-method argument lists
#' @param with_value Constant value for imp_method="with"
#' @param comparisons Comparisons for DE. If NULL, generates all pairwise
#' @param control Control condition. If NULL, compares all
#' @param logFC_threshold LogFC threshold for significance (default: 0)
#' @param alpha Adjusted p-value threshold (default: 0.05)
#' @param eBayes_trend Use trend estimation in eBayes. If NULL (default), it is
#'   resolved from de_method: TRUE for "limma", FALSE for "limpa".
#' @param eBayes_robust Use robust estimation in eBayes. If NULL (default), it is
#'   resolved from de_method: TRUE for "limma", FALSE for "limpa".
#' @param de_method DE method: "limma" (default) or "limpa" (probabilistic, requires imp_method="limpa")
#' @param covariate_df Data frame with Column + covariate column(s) for paired/blocked design (default: NULL)
#' @param covariate_column Name(s) of the covariate column(s) for the DE model.
#'   Single string (e.g., "Subject") or character vector (e.g., c("Subject", "Batch")). Default: NULL
#' @param bio_replicate_column Column name in covariate_df identifying biological replicates
#'   (e.g., "Patient", "Subject"). Used with limma::duplicateCorrelation() to account for
#'   technical replicates or paired designs via random effect blocking. Default: NULL
#' @param export_normalized Export normalized matrix (default: TRUE)
#' @param export_imputed Export imputed matrix (default: TRUE)
#' @param export_format Export format: "tsv", "parquet", or "both" (default: "tsv")
#' @param export_volcano Export VolcanoPlot_Input (default: TRUE)
#' @param export_boxplot Export BoxPlot_Input (default: TRUE)
#' @param export_pca Export PCA_Input (default: TRUE)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return List of class proteomics_result with:
#'   \itemize{
#'     \item se_proc: Processed SummarizedExperiment with all assays
#'     \item DEPs_results: Data frame with differential expression results
#'     \item BoxPlot_Input: Long-format data frame with one row per protein,
#'       sample and assay (columns Column, Assay, Intensity, Condition,
#'       Replicate, Protein.IDs). Feed it to [boxplot_highchart_list()].
#'     \item PCA_Input: Long-format data frame with the imputed intensities plus
#'       one adjP_* column per comparison and a sig_any flag. Feed it to
#'       [pca_highchart_list()], [proteomics_heatmap()] or
#'       [proteomics_heatmap_list()].
#'     \item comparisons: Comparisons performed
#'     \item parameters: Parameters used
#'   }
#'
#'   `BoxPlot_Input` and `PCA_Input` are always returned, whether or not
#'   `export_dir` is given: the plotting layer needs them, and requiring a round
#'   trip through disk to obtain them would be gratuitous.
#'
#' @examples
#' data(nadia_dia)
#'
#' # Normalization, imputation and differential expression in one call
#' res <- process_proteomics(nadia_dia, verbose = FALSE)
#' res
#'
#' SummarizedExperiment::assayNames(res$se_proc)
#' table(res$DEPs_results$Comparison, res$DEPs_results$Change)
#'
#' # The plotting layer is fed straight from the result, no files involved
#' names(res$BoxPlot_Input)
#' names(res$PCA_Input)
#'
#' # Another combination: quantile normalization and a single imputation method
#' res2 <- process_proteomics(nadia_dia, norm_method = "quantile",
#'                            imp_method = "QRILC", verbose = FALSE)
#' SummarizedExperiment::assayNames(res2$se_proc)
#'
#' @export
process_proteomics <- function(
    preprocessing,
    export_dir = NULL,
    min_reps_filter = NULL,
    min_groups_filter = 1,
    norm_method = "cycloess",
    cyclic_loess_method = "fast",
    cyclic_loess_iterations = 3,
    cyclic_loess_span = 0.7,
    batch_correct = FALSE,
    batch_column = "Batch",
    batch_algorithm = "ComBat",
    batch_ComBat_mode = 1,
    batch_covariates = NULL,
    batch_qualitycontrol = FALSE,
    imp_method = "combo",
    mar_method = "Impseqrob",
    mnar_method = "min",
    prop_na_mnar = 0.51,
    prop_present_mar = 0.5,
    min_present_mar = 1,
    require_n_conditions = 1,
    max_na_prop = NULL,
    method_args = list(),
    with_value = NA_real_,
    comparisons = NULL,
    control = NULL,
    logFC_threshold = 0,
    alpha = 0.05,
    eBayes_trend = NULL,
    eBayes_robust = NULL,
    de_method = "limma",
    covariate_df = NULL,
    covariate_column = NULL,
    bio_replicate_column = NULL,
    export_normalized = TRUE,
    export_imputed = TRUE,
    export_format = "tsv",
    export_volcano = TRUE,
    export_boxplot = TRUE,
    export_pca = TRUE,
    verbose = TRUE
) {
  # =========================================================================
  # VALIDATIONS
  # =========================================================================

  if (!inherits(preprocessing, "proteomics_data")) {
    stop("The 'preprocessing' argument must be the result of preprocess_spectronaut() or preprocess_tmt()")
  }

  # Without export_dir nothing is written to disk. The function must not create
  # files or directories in the user's workspace unless explicitly told where
  # (a Bioconductor requirement). Turning the flags off here is enough to cover
  # all the later export blocks.
  if (is.null(export_dir)) {
    export_normalized <- FALSE
    export_imputed    <- FALSE
    export_volcano    <- FALSE
    export_boxplot    <- FALSE
    export_pca        <- FALSE
  } else if (!dir.exists(export_dir)) {
    dir.create(export_dir, recursive = TRUE)
  }

  # =========================================================================
  # 1. PREPARE DATA FROM PREPROCESSING
  # =========================================================================

  if (verbose) message("=== PREPARING DATA ===")

  metadata <- .prepare_metadata(preprocessing, covariate_df = covariate_df)
  protein_data <- .prepare_protein_data(preprocessing)

  if (verbose) {
    message("- Metadata: ", nrow(metadata), " samples")
    message("- Proteins: ", nrow(protein_data), " initial proteins")
  }

  # =========================================================================
  # 2. NORMALIZATION (Normalization.R)
  # =========================================================================

  norm_result <- normalize_proteomics(
    data = protein_data,
    metadata = metadata,
    min_reps = min_reps_filter,
    min_groups = min_groups_filter,
    norm_method = norm_method,
    cyclic_loess_method = cyclic_loess_method,
    cyclic_loess_iterations = cyclic_loess_iterations,
    cyclic_loess_span = cyclic_loess_span,
    verbose = verbose
  )

  se <- norm_result$se

  # Export normalized matrix (skip when norm_method = "log2": no extra assay)
  if (export_normalized && norm_method != "log2") {
    x_norm <- SummarizedExperiment::assay(se, norm_method)
    norm_file <- file.path(export_dir, paste0("matrix_log2_", norm_method, ".tsv"))
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_tsv(
        data.frame(ProteinGroups = rownames(x_norm), x_norm, check.names = FALSE),
        norm_file
      )
    } else {
      write.table(
        data.frame(ProteinGroups = rownames(x_norm), x_norm, check.names = FALSE),
        norm_file, sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
    if (verbose) message("- Exported: ", basename(norm_file))
  }

  # =========================================================================
  # 2b. BATCH CORRECTION (optional -- Batch_Correction.R / BERT)
  # =========================================================================

  input_to_imputation <- norm_method

  if (batch_correct) {
    if (!batch_column %in% colnames(SummarizedExperiment::colData(se))) {
      stop("batch_correct=TRUE but batch column '", batch_column,
           "' not found in colData(se).\n",
           "  Ensure covariate_df contains a '", batch_column, "' column.")
    }

    se <- batch_correct_proteomics(
      se                   = se,
      assay_name           = norm_method,
      batch_column         = batch_column,
      corrected_assay_name = "BERT",
      algorithm            = batch_algorithm,
      ComBat_mode          = batch_ComBat_mode,
      covariates           = batch_covariates,
      qualitycontrol       = batch_qualitycontrol,
      verbose              = verbose
    )

    input_to_imputation <- "BERT"

    # Export batch-corrected matrix
    if (export_normalized) {
      x_bc <- SummarizedExperiment::assay(se, "BERT")
      bc_file <- file.path(export_dir,
                           paste0("matrix_log2_", norm_method, "_BERT.tsv"))
      if (requireNamespace("readr", quietly = TRUE)) {
        readr::write_tsv(
          data.frame(ProteinGroups = rownames(x_bc), x_bc, check.names = FALSE),
          bc_file
        )
      } else {
        write.table(
          data.frame(ProteinGroups = rownames(x_bc), x_bc, check.names = FALSE),
          bc_file, sep = "\t", quote = FALSE, row.names = FALSE
        )
      }
      if (verbose) message("- Exported: ", basename(bc_file))
    }
  }

  # =========================================================================
  # 3. IMPUTATION (Imputation.R)
  # =========================================================================

  imp_result <- impute_proteomics(
    se = se,
    normalized_assay_name = input_to_imputation,
    imputed_assay_name = NULL,
    imp_method = imp_method,
    mar_method = mar_method,
    mnar_method = mnar_method,
    prop_na_mnar = prop_na_mnar,
    prop_present_mar = prop_present_mar,
    min_present_mar = min_present_mar,
    require_n_conditions = require_n_conditions,
    max_na_prop = max_na_prop,
    method_args = method_args,
    with_value = with_value,
    verbose = verbose
  )

  # Derive assay_label from the imputed SE
  assay_label <- setdiff(
    SummarizedExperiment::assayNames(imp_result$se),
    SummarizedExperiment::assayNames(se)
  )
  if (length(assay_label) == 0) {
    # imp_method="none" adds no new assay; fall back to normalized
    assay_label <- norm_method
  } else {
    assay_label <- assay_label[1]
  }

  se_proc <- imp_result$se

  # Export imputed matrix (use raw matrix, matching pre-split behavior)
  if (export_imputed) {
    x_imputed_export <- imp_result$x_imputed
    imp_file <- file.path(export_dir, paste0("matrix_log2_", norm_method, "_", assay_label, ".tsv"))
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_tsv(
        data.frame(ProteinGroups = rownames(x_imputed_export),
                   x_imputed_export, check.names = FALSE),
        imp_file
      )
    } else {
      write.table(
        data.frame(ProteinGroups = rownames(x_imputed_export),
                   x_imputed_export, check.names = FALSE),
        imp_file, sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
    if (verbose) message("- Exported: ", basename(imp_file))
  }

  # =========================================================================
  # 4. DIFFERENTIAL EXPRESSION ANALYSIS (DEAnalysis.R)
  # =========================================================================

  de_result <- de_analysis_proteomics(
    se = se_proc,
    assay_name = assay_label,
    comparisons = comparisons,
    control = control,
    logFC_threshold = logFC_threshold,
    alpha = alpha,
    p_adj = TRUE,
    eBayes_trend = eBayes_trend,
    eBayes_robust = eBayes_robust,
    de_method = de_method,
    covariate_column = covariate_column,
    bio_replicate_column = bio_replicate_column,
    condition_column = "Condition",
    verbose = verbose
  )

  DEPs_results <- de_result$DEPs_results
  comparisons <- de_result$comparisons

  # 4b. Add MissGlobal, MissComp, MissCND1, MissCND2
  DEPs_results <- .compute_missing_pct(
    se = se_proc,
    DEPs_results = DEPs_results,
    norm_assay_name = norm_method
  )

  # =========================================================================
  # 5. BUILD THE VISUALIZATION INPUTS (and export them if asked to)
  # =========================================================================

  # These two tables are what the plotting layer consumes. They are built
  # unconditionally and returned, so that a caller can pipe straight into
  # boxplot_highchart_list(), pca_highchart_list() or proteomics_heatmap()
  # without having to write anything to disk. Building them is cheap (well under
  # a tenth of a second on a typical dataset).
  boxplot_data <- .se_to_long(se_proc, assay_names = c("log2", assay_label))
  pca_data     <- .prepare_pca_input(se_proc, DEPs_results, assay_label, alpha)

  if (export_volcano || export_boxplot || export_pca) {
    if (verbose) message("\n=== EXPORTING FILES FOR VISUALIZATION ===")

    # Suffix with norm_method and assay_label (avoid redundancy if they are equal)
    viz_suffix <- if (identical(assay_label, norm_method)) {
      paste0("_", norm_method)
    } else {
      paste0("_", norm_method, "_", assay_label)
    }

    # 5.1 VolcanoPlot_Input (DEPs_results)
    if (export_volcano) {
      volcano_file <- file.path(export_dir, paste0("VolcanoPlot_Input", viz_suffix))
      .export_data(DEPs_results, volcano_file, export_format)
      if (verbose) message("- VolcanoPlot_Input exported")
    }

    # 5.2 BoxPlot_Input (SE -> long format)
    if (export_boxplot) {
      boxplot_file <- file.path(export_dir, paste0("BoxPlot_Input", viz_suffix))
      .export_data(boxplot_data, boxplot_file, export_format)
      if (verbose) message("- BoxPlot_Input exported")
    }

    # 5.3 PCA_Input (SE + DE in long format)
    if (export_pca) {
      pca_file <- file.path(export_dir, paste0("PCA_Input", viz_suffix))
      .export_data(pca_data, pca_file, export_format)
      if (verbose) message("- PCA_Input exported")
    }
  }

  # =========================================================================
  # 6. RETURN RESULT
  # =========================================================================

  result <- list(
    se_proc = se_proc,
    DEPs_results = DEPs_results,
    BoxPlot_Input = boxplot_data,
    PCA_Input = pca_data,
    comparisons = comparisons,
    parameters = list(
      min_reps_filter = norm_result$filter_summary$min_reps,
      min_groups_filter = min_groups_filter,
      norm_method = norm_method,
      cyclic_loess_method = cyclic_loess_method,
      cyclic_loess_iterations = cyclic_loess_iterations,
      cyclic_loess_span = cyclic_loess_span,
      batch_correct = batch_correct,
      batch_column = batch_column,
      batch_algorithm = batch_algorithm,
      batch_ComBat_mode = batch_ComBat_mode,
      batch_covariates = batch_covariates,
      batch_qualitycontrol = batch_qualitycontrol,
      imp_method = imp_method,
      mar_method = mar_method,
      mnar_method = mnar_method,
      prop_na_mnar = prop_na_mnar,
      prop_present_mar = prop_present_mar,
      max_na_prop = max_na_prop,
      logFC_threshold = logFC_threshold,
      alpha = alpha,
      eBayes_trend = eBayes_trend,
      eBayes_robust = eBayes_robust,
      de_method = de_method,
      covariate_column = covariate_column,
      bio_replicate_column = bio_replicate_column,
      export_dir = export_dir,
      export_format = export_format,
      export_volcano = export_volcano,
      export_boxplot = export_boxplot,
      export_pca = export_pca
    )
  )

  class(result) <- c("proteomics_result", "list")

  if (verbose) message("\n=== PROCESSING COMPLETED ===")

  result
}

# =============================================================================
# PRINT METHOD
# =============================================================================

#' Print method for proteomics_result
#'
#' @param x proteomics_result object
#' @param ... Additional arguments (ignored)
#' @return `x`, invisibly. Called for the summary it prints: the assays of the
#'   SummarizedExperiment, the comparisons performed, the number of differential
#'   proteins and the parameters used.
#'
#' @examples
#' data(nadia_dia)
#' res <- process_proteomics(nadia_dia, verbose = FALSE)
#' print(res)
#'
#' @export
print.proteomics_result <- function(x, ...) {
  cat("=== Proteomics Processing Result ===\n\n")

  # SE summary
  se <- x$se_proc
  cat("SummarizedExperiment:\n")
  cat("  - Proteins:", nrow(se), "\n")
  cat("  - Samples:", ncol(se), "\n")
  cat("  - Assays:", paste(SummarizedExperiment::assayNames(se), collapse = ", "), "\n")

  # Condition summary
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if ("Condition" %in% names(cd)) {
    cat("  - Conditions:", paste(unique(cd$Condition), collapse = ", "), "\n")
  }

  cat("\nDifferential Results:\n")
  cat("  - Total rows:", nrow(x$DEPs_results), "\n")
  cat("  - Comparisons:", paste(unique(x$DEPs_results$Comparison), collapse = ", "), "\n")

  for (comp in unique(x$DEPs_results$Comparison)) {
    subset <- x$DEPs_results[x$DEPs_results$Comparison == comp, ]
    n_up <- sum(subset$Change == "Up", na.rm = TRUE)
    n_down <- sum(subset$Change == "Down", na.rm = TRUE)
    cat("    ", comp, ": Up=", n_up, ", Down=", n_down, "\n", sep = "")
  }

  cat("\nParameters:\n")
  cat("  - Normalization:", x$parameters$norm_method, "\n")
  if (identical(x$parameters$norm_method, "cycloess")) {
    cat("    - Cyclic Loess method:", x$parameters$cyclic_loess_method, "\n")
    cat("    - Cyclic Loess iterations:", x$parameters$cyclic_loess_iterations, "\n")
    cat("    - Cyclic Loess span:", x$parameters$cyclic_loess_span, "\n")
  }
  cat("  - Imputation:", x$parameters$imp_method, "\n")
  if (identical(x$parameters$imp_method, "combo")) {
    cat("    - MAR method:", x$parameters$mar_method, "\n")
    cat("    - MNAR method:", x$parameters$mnar_method, "\n")
  }
  cat("  - Alpha:", x$parameters$alpha, "\n")
  cat("  - logFC threshold:", x$parameters$logFC_threshold, "\n")
  cat("  - Output directory:",
      x$parameters$export_dir %||% "(no export)", "\n")

  invisible(x)
}
