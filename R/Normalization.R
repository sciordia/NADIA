# =============================================================================
# Normalization Module
# =============================================================================
#
# Functions for proteomics data normalization:
#   - Zero-to-NA conversion
#   - Protein filtering by group presence
#   - SummarizedExperiment creation
#   - 13 normalization methods (cycloess default)
#
# Dependencies (required):
#   - SummarizedExperiment, S4Vectors
#
# Dependencies (optional, per method):
#   - limma : cycloess, quantile
#   - MASS  : Rlr
#   - vsn   : vsn
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================

# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Convert zero values to NA
#'
#' Replaces 0 values in a matrix or data frame with NA.
#' Useful for DIA-NN data where 0 indicates non-detection.
#'
#' @param data Numeric matrix or data frame
#' @return Object of same type with 0 converted to NA
#' @keywords internal
.zero_to_missing <- function(data) {
  if (is.data.frame(data)) {
    numeric_cols <- vapply(data, is.numeric, logical(1))
    data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
      x[x == 0] <- NA
      x
    })
  } else if (is.matrix(data)) {
    data[data == 0] <- NA
  }
  data
}

#' Filter proteins by group presence
#'
#' Keeps proteins that have enough non-NA values in at least
#' a minimum number of experimental groups.
#'
#' @param data Intensity matrix (proteins x samples)
#' @param metadata Data frame with sample information
#' @param min_reps Minimum replicates with non-NA values per group.
#'   If NULL, uses half of the smallest group size.
#' @param min_groups Minimum groups meeting min_reps (default: 1)
#' @param grouping_column Name of grouping column in metadata (default: "Condition")
#' @return List with:
#'   - data: Filtered matrix
#'   - keep: Logical vector of kept rows
#'   - summary: Filtering summary
#' @keywords internal
.filter_proteins_by_group <- function(
    data,
    metadata,
    min_reps = NULL,
    min_groups = 1,
    grouping_column = "Condition"
) {
  # Validations
  stopifnot(is.matrix(data) || is.data.frame(data))
  data <- as.matrix(data)

  if (!grouping_column %in% names(metadata)) {
    stop("Column '", grouping_column, "' does not exist in metadata")
  }

  # Align samples
  if (!is.null(rownames(metadata))) {
    common_samples <- intersect(colnames(data), rownames(metadata))
    if (length(common_samples) == 0) {
      stop("No samples in common between data and metadata")
    }
    data <- data[, common_samples, drop = FALSE]
    metadata <- metadata[common_samples, , drop = FALSE]
  }

  groups <- as.factor(metadata[[grouping_column]])
  group_levels <- levels(groups)

  # Auto-compute min_reps if not specified
  if (is.null(min_reps)) {
    group_sizes <- table(groups)
    min_reps <- max(1, floor(min(group_sizes) / 2))
  }

  # Compute non-NA counts per group for each protein
  n_present_per_group <- vapply(group_levels, function(g) {
    cols <- which(groups == g)
    if (length(cols) == 0) return(rep(0, nrow(data)))
    rowSums(!is.na(data[, cols, drop = FALSE]))
  }, numeric(nrow(data)))

  if (!is.matrix(n_present_per_group)) {
    n_present_per_group <- matrix(n_present_per_group, ncol = 1)
  }

  # Count groups meeting criteria
  groups_ok <- rowSums(n_present_per_group >= min_reps)
  keep <- groups_ok >= min_groups

  list(
    data = data[keep, , drop = FALSE],
    keep = keep,
    n_present_per_group = n_present_per_group,
    summary = list(
      n_total = nrow(data),
      n_keep = sum(keep),
      n_drop = sum(!keep),
      min_reps = min_reps,
      min_groups = min_groups
    )
  )
}

#' Create SummarizedExperiment from proteomics data
#'
#' @param data Data frame with proteins and intensity values
#' @param metadata Data frame with sample information
#' @param protein_column Protein ID column name
#' @param gene_column Gene name column name
#' @param condition_column Condition column name in metadata
#' @param label_column Sample label column name in metadata
#' @return SummarizedExperiment with assays: raw, log2
#' @keywords internal
.load_proteomics_data <- function(
    data,
    metadata,
    protein_column = "ProteinGroups",
    gene_column = "GeneNames",
    condition_column = "Condition",
    label_column = "Column"
) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required")
  }

  # Validate required columns
  if (!protein_column %in% names(data)) {
    stop("Column '", protein_column, "' not found in data")
  }
  if (!label_column %in% names(metadata)) {
    stop("Column '", label_column, "' not found in metadata")
  }
  if (!condition_column %in% names(metadata)) {
    stop("Column '", condition_column, "' not found in metadata")
  }

  # Identify intensity columns
  annotation_cols <- c(protein_column, gene_column, "UniqPepts")
  annotation_cols <- intersect(annotation_cols, names(data))
  intensity_cols <- setdiff(names(data), annotation_cols)

  # Extract matrices
  annotation <- data[, annotation_cols, drop = FALSE]
  intensity <- as.matrix(data[, intensity_cols, drop = FALSE])
  storage.mode(intensity) <- "double"

  # Unique IDs
  protein_ids <- as.character(annotation[[protein_column]])
  protein_ids <- make.unique(protein_ids)
  rownames(intensity) <- protein_ids
  rownames(annotation) <- protein_ids

  # Align metadata with intensity columns
  rownames(metadata) <- metadata[[label_column]]
  common_samples <- intersect(colnames(intensity), rownames(metadata))

  if (length(common_samples) == 0) {
    stop("No samples in common between data and metadata")
  }

  intensity <- intensity[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  # Create rowData
  row_data <- S4Vectors::DataFrame(annotation)
  names(row_data)[names(row_data) == protein_column] <- "Protein.IDs"
  if (gene_column %in% names(row_data)) {
    names(row_data)[names(row_data) == gene_column] <- "Gene.Names"
  }
  row_data$IDs <- rownames(row_data)

  # Create colData
  col_data <- S4Vectors::DataFrame(metadata)

  # Create assays: raw and log2
  raw_assay <- intensity
  log2_assay <- log2(intensity)
  log2_assay[is.infinite(log2_assay)] <- NA

  # Create SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = raw_assay, log2 = log2_assay),
    rowData = row_data,
    colData = col_data,
    metadata = list(
      condition = condition_column,
      label = label_column
    )
  )

  se
}

#' NA overview for a SummarizedExperiment
#'
#' @param se SummarizedExperiment
#' @param assay_name Assay name to evaluate (default: "log2")
#' @return Data frame with NA statistics per sample
#' @keywords internal
.get_NA_overview <- function(se, assay_name = "log2") {
  stopifnot(inherits(se, "SummarizedExperiment"))

  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", assay_name, "' not found")
  }

  x <- SummarizedExperiment::assay(se, assay_name)

  # Per sample
  total_vals <- nrow(x)
  na_counts <- colSums(is.na(x))
  na_pct <- na_counts / total_vals * 100

  sample_stats <- data.frame(
    Sample = colnames(x),
    Total.Values = total_vals,
    NA.Values = na_counts,
    NA.Percentage = round(na_pct, 2),
    stringsAsFactors = FALSE
  )

  # Global
  total_cells <- length(x)
  total_na <- sum(is.na(x))

  attr(sample_stats, "global") <- list(
    Total.Cells = total_cells,
    Total.NA = total_na,
    NA.Percentage = round(total_na / total_cells * 100, 2)
  )

  sample_stats
}

# =============================================================================
# INTERNAL NORMALIZATION FUNCTIONS
# =============================================================================
#
# All functions return a numeric matrix in log2 scale with the same
# rownames/colnames as the input.
#
# Group A (.norm_log2norm .. .norm_vsn): receive x_raw (linear).
# Group B (.norm_quantile .. .norm_quantile_robust): receive x_log2.
#
# Note: eqmedians belongs to Group A (receives x_raw) but applies log2() internally.

# --- Group A: x_raw -> log2 ---

.norm_log2norm <- function(x_raw) {
  x <- log2(x_raw)
  x[is.infinite(x)] <- NA
  x
}

# NOTE: "GlobalMedian" (.norm_ginorm) and "GlobalMean" (.norm_globalmean)
# normalize by the column SUM, rescaled to the median / mean of those sums.
# They do not use the column median/mean (that is medianNorm/meanNorm). Both
# methods differ only in the global constant log2(median(S)) vs log2(mean(S));
# since every downstream metric (variance, correlation, centred PCA/MDS,
# PCV/PMAD/PEV) is invariant to a global shift, they produce practically
# identical results in the benchmark. Both are kept for backward
# compatibility, but the redundancy is documented here.
.norm_ginorm <- function(x_raw) {
  col_sums <- colSums(x_raw, na.rm = TRUE)
  x <- log2(sweep(x_raw, 2, col_sums / median(col_sums), "/"))
  x[is.infinite(x)] <- NA
  x
}

.norm_globalmean <- function(x_raw) {
  col_sums <- colSums(x_raw, na.rm = TRUE)
  x <- log2(sweep(x_raw, 2, col_sums / mean(col_sums), "/"))
  x[is.infinite(x)] <- NA
  x
}

# These two receive x_raw but apply log2 internally before centering.

.norm_eqmedians <- function(x_raw) {
  x <- log2(x_raw)
  x[is.infinite(x)] <- NA
  col_medians <- apply(x, 2, median, na.rm = TRUE)
  # MSstats equalizeMedians: x - colMedian + median(colMedians)
  sweep(x, 2, col_medians - median(col_medians), "-")
}


.norm_vsn <- function(x_raw) {
  if (!requireNamespace("vsn", quietly = TRUE)) {
    stop("Package 'vsn' is required for the vsn method. ",
         "Install it with BiocManager::install('vsn')")
  }
  vsn::justvsn(x_raw)
}

# --- Group B: x_log2 -> log2 ---

.norm_quantile <- function(x_log2) {
  # limma::normalizeQuantiles handles NA consistently (it interpolates the
  # per-column quantiles onto a common grid), unlike
  # preprocessCore::normalize.quantiles, which propagates NaN with DIA data.
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for the quantile method. ",
         "Install it with BiocManager::install('limma')")
  }
  res <- limma::normalizeQuantiles(x_log2)
  dimnames(res) <- dimnames(x_log2)
  res
}

.norm_rlr <- function(x_log2) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for the Rlr method. ",
         "Install it with install.packages('MASS')")
  }
  row_medians <- apply(x_log2, 1, median, na.rm = TRUE)
  x <- x_log2
  for (j in seq_len(ncol(x_log2))) {
    col   <- x_log2[, j]
    valid <- is.finite(col) & is.finite(row_medians)
    if (sum(valid) < 2L) next
    fit <- tryCatch(MASS::rlm(col[valid] ~ row_medians[valid]),
                    error = function(e) NULL)
    if (is.null(fit)) next
    intercept <- coef(fit)[1L]
    slope     <- coef(fit)[2L]
    if (is.na(slope) || abs(slope) < .Machine$double.eps) next
    x[, j] <- (col - intercept) / slope
  }
  x
}

.norm_mad <- function(x_log2) {
  grand_median <- median(x_log2, na.rm = TRUE)
  grand_mad    <- mad(x_log2,    na.rm = TRUE)
  col_medians  <- apply(x_log2, 2, median, na.rm = TRUE)
  col_mads     <- apply(x_log2, 2, mad,    na.rm = TRUE)
  col_mads[col_mads == 0] <- 1
  x <- sweep(x_log2, 2, col_medians, "-")
  x <- sweep(x,      2, col_mads,    "/")
  x * grand_mad + grand_median
}

.norm_mediannorm <- function(x_raw) {
  # NormalyzerDE/PRONE: (x_raw / colMedian) * mean(colMedians) -> log2
  col_medians <- apply(x_raw, 2, median, na.rm = TRUE)
  x <- log2(sweep(x_raw, 2, col_medians / mean(col_medians), "/"))
  x[is.infinite(x)] <- NA
  x
}

.norm_meannorm <- function(x_raw) {
  # NormalyzerDE/PRONE: (x_raw / colMean) * mean(colMeans) -> log2
  col_means <- colMeans(x_raw, na.rm = TRUE)
  x <- log2(sweep(x_raw, 2, col_means / mean(col_means), "/"))
  x[is.infinite(x)] <- NA
  x
}

.norm_quantile_robust <- function(x_log2) {
  # Robust quantile normalization -- base R, no external dependencies.
  # Uses the median (instead of the mean) of the sorted values as the
  # reference distribution, which makes it robust to samples with extreme
  # values.
  #
  # With NA, each column is first interpolated onto a common grid of n_row
  # quantiles in [0,1] before taking the per-row median; this prevents mixing
  # different quantiles across columns with a different number of observations
  # (the sort(na.last=TRUE) bug: position i was quantile i/m, not i/n).
  n_row <- nrow(x_log2)
  n_col <- ncol(x_log2)
  x_norm <- x_log2
  grid   <- if (n_row > 1L) (seq_len(n_row) - 1L) / (n_row - 1L) else 0

  # Reference: median per quantile position across the interpolated columns.
  interp_cols <- vapply(seq_len(n_col), function(j) {
    obs <- sort(x_log2[!is.na(x_log2[, j]), j])
    m   <- length(obs)
    if (m == 0L) return(rep(NA_real_, n_row))
    if (m == 1L) return(rep(obs, n_row))
    approx((seq_len(m) - 1L) / (m - 1L), obs, xout = grid, rule = 2L)$y
  }, numeric(n_row))
  ref <- apply(interp_cols, 1, median, na.rm = TRUE)

  # Per-column mapping: value rank -> grid position -> reference value.
  for (j in seq_len(n_col)) {
    col   <- x_log2[, j]
    valid <- !is.na(col)
    n_j   <- sum(valid)
    if (n_j == 0L) next
    r       <- rank(col[valid], ties.method = "average")
    ref_pos <- if (n_j > 1L) (r - 1) / (n_j - 1L) else 0
    x_norm[valid, j] <- approx(grid, ref, xout = ref_pos, rule = 2L)$y
  }
  x_norm
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Normalize proteomics data
#'
#' Complete normalization pipeline: zero-to-NA conversion, protein filtering
#' by group presence, SummarizedExperiment creation, and normalization with
#' the selected method. Output is always in log2 scale.
#'
#' @param data Either a `proteomics_data` object -- the output of
#'   [preprocess_spectronaut()], [preprocess_tmt()] or [preprocess_lfq()] -- in
#'   which case `metadata` is derived from it and must be left `NULL`, or a data
#'   frame with ProteinGroups, GeneNames, UniqPepts and the intensity columns.
#' @param metadata Data frame with Column, Condition, Replicate. Only needed when
#'   `data` is a data frame.
#' @param min_reps Minimum replicates for filtering. NULL = auto: floor(min_group_size / 2)
#' @param min_groups Minimum groups meeting min_reps (default: 1)
#' @param norm_method Normalization method (default: "cycloess"). One of:
#'   \itemize{
#'     \item Group A (input: raw intensities): "log2Norm", "GlobalMedian",
#'       "GlobalMean", "eqmedians", "vsn"
#'     \item Group B (input: log2 assay): "log2" (no extra normalization),
#'       "quantile", "Rlr", "MAD", "cycloess",
#'       "medianNorm", "meanNorm", "quantile.robust"
#'   }
#' @param cyclic_loess_method Cyclic Loess method: "fast" or "pairs" (default: "fast")
#' @param cyclic_loess_iterations Number of iterations for Cyclic Loess (default: 3)
#' @param cyclic_loess_span Span parameter for Cyclic Loess (default: 0.7)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return List with:
#'   \itemize{
#'     \item se: SummarizedExperiment with assays raw, log2, and <norm_method>
#'       (assays raw + log2 only when norm_method = "log2")
#'     \item filter_summary: Filtering summary
#'     \item na_overview: NA statistics
#'   }
#'
#' @examples
#' \dontrun{
#' norm_result <- normalize_proteomics(
#'   data = protein_data,
#'   metadata = metadata,
#'   norm_method = "cycloess",
#'   cyclic_loess_method = "fast"
#' )
#' norm_result2 <- normalize_proteomics(
#'   data = protein_data,
#'   metadata = metadata,
#'   norm_method = "GlobalMedian"
#' )
#' }
#'
#' @export
normalize_proteomics <- function(
    data,
    metadata                = NULL,
    min_reps                = NULL,
    min_groups              = 1,
    norm_method             = "cycloess",
    cyclic_loess_method     = c("fast", "pairs"),
    cyclic_loess_iterations = 3,
    cyclic_loess_span       = 0.7,
    verbose                 = TRUE
) {
  # Accept the object returned by preprocess_spectronaut()/_tmt()/_lfq() directly.
  # Otherwise the function would be unusable on its own, since the two tables it
  # expects are only produced by internal helpers.
  if (inherits(data, "proteomics_data")) {
    if (!is.null(metadata)) {
      stop("When 'data' is a proteomics_data object, 'metadata' must not be ",
           "supplied: it is taken from the object itself.")
    }
    preprocessing <- data
    metadata      <- .prepare_metadata(preprocessing)
    data          <- .prepare_protein_data(preprocessing)
  } else if (is.null(metadata)) {
    stop("'metadata' is required when 'data' is a data frame. Pass the object ",
         "returned by preprocess_spectronaut(), preprocess_tmt() or ",
         "preprocess_lfq() to have both derived automatically.")
  }

  # Match cycloess sub-parameters
  cyclic_loess_method <- match.arg(cyclic_loess_method)

  # Validate required packages
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required. ",
         "Install it with BiocManager::install('SummarizedExperiment')")
  }

  # =========================================================================
  # 1. ZERO TO NA CONVERSION
  # =========================================================================

  if (verbose) cat("\n=== CONVERTING ZEROS TO NA ===\n")

  annotation_cols <- c("ProteinGroups", "GeneNames", "UniqPepts")
  intensity_cols <- setdiff(names(data), annotation_cols)

  intensity_mat <- as.matrix(data[, intensity_cols])
  n_zeros <- sum(intensity_mat == 0, na.rm = TRUE)
  intensity_mat <- .zero_to_missing(intensity_mat)
  rownames(intensity_mat) <- data$ProteinGroups

  if (verbose) cat("- Zeros converted to NA:", n_zeros, "\n")

  # =========================================================================
  # 2. FILTER PROTEINS BY GROUP PRESENCE
  # =========================================================================

  if (verbose) cat("\n=== FILTERING PROTEINS BY PRESENCE ===\n")

  filtered <- .filter_proteins_by_group(
    data = intensity_mat,
    metadata = metadata,
    min_reps = min_reps,
    min_groups = min_groups,
    grouping_column = "Condition"
  )

  if (verbose) {
    cat("- Proteins before:", filtered$summary$n_total, "\n")
    cat("- Proteins after:", filtered$summary$n_keep, "\n")
    cat("- Proteins removed:", filtered$summary$n_drop, "\n")
    cat("- Min replicates:", filtered$summary$min_reps, "\n")
    cat("- Min groups:", filtered$summary$min_groups, "\n")
  }

  # Filter protein_data
  protein_data_filtered <- data[filtered$keep, , drop = FALSE]
  # Propagate zero-to-NA conversion into the matrix that feeds the SE
  # (intensity_mat already has zeros -> NA; aligned by position with filtered$keep)
  protein_data_filtered[, intensity_cols] <- intensity_mat[filtered$keep, , drop = FALSE]

  # =========================================================================
  # 3. CREATE SUMMARIZEDEXPERIMENT
  # =========================================================================

  if (verbose) cat("\n=== CREATING SUMMARIZEDEXPERIMENT ===\n")

  se <- .load_proteomics_data(
    data = protein_data_filtered,
    metadata = metadata,
    protein_column = "ProteinGroups",
    gene_column = "GeneNames",
    condition_column = "Condition",
    label_column = "Column"
  )

  na_overview <- NULL
  if (verbose) {
    na_overview <- .get_NA_overview(se, "log2")
    global_na <- attr(na_overview, "global")
    cat("- Global NA:", global_na$NA.Percentage, "%\n")
  }

  # =========================================================================
  # 4. NORMALIZATION
  # =========================================================================

  # Methods that accept x_raw as input.
  # (eqmedians applies log2 internally but still receives x_raw to stay
  #  consistent with the other Group A methods.)
  .raw_methods <- c(
    "log2Norm", "GlobalMedian", "GlobalMean",
    "eqmedians",
    "vsn",
    "medianNorm", "meanNorm"
  )

  norm_method <- match.arg(norm_method, c(
    "log2Norm", "eqmedians", "GlobalMedian", "GlobalMean",
    "log2", "quantile", "Rlr", "MAD", "cycloess",
    "medianNorm", "meanNorm",
    "vsn", "quantile.robust"
  ))

  x_raw  <- SummarizedExperiment::assay(se, "raw")
  x_log2 <- SummarizedExperiment::assay(se, "log2")

  if (verbose) cat("\n=== NORMALIZING (method:", norm_method, ") ===\n")

  if (norm_method == "log2") {
    if (verbose) cat("- Method 'log2': no additional normalization\n")
  } else {
    x_input <- if (norm_method %in% .raw_methods) x_raw else x_log2

    x_norm <- switch(norm_method,
      "log2Norm"        = .norm_log2norm(x_input),
      "GlobalMedian"    = .norm_ginorm(x_input),
      "GlobalMean"      = .norm_globalmean(x_input),
      "eqmedians"       = .norm_eqmedians(x_input),
      "quantile"        = .norm_quantile(x_input),
      "Rlr"             = .norm_rlr(x_input),
      "MAD"             = .norm_mad(x_input),
      "cycloess"        = {
                            if (!requireNamespace("limma", quietly = TRUE))
                              stop("Package 'limma' is required for the cycloess method. ",
                                   "Install it with BiocManager::install('limma')")
                            limma::normalizeCyclicLoess(
                              x_input,
                              method     = cyclic_loess_method,
                              iterations = cyclic_loess_iterations,
                              span       = cyclic_loess_span)
                          },
      "medianNorm"        = .norm_mediannorm(x_input),
      "meanNorm"          = .norm_meannorm(x_input),
      "vsn"               = .norm_vsn(x_input),
      "quantile.robust"   = .norm_quantile_robust(x_input)
    )

    rownames(x_norm) <- rownames(x_log2)
    colnames(x_norm) <- colnames(x_log2)
    SummarizedExperiment::assay(se, norm_method) <- x_norm
  }

  if (verbose) cat("- Available assays:",
                   paste(SummarizedExperiment::assayNames(se), collapse = ", "), "\n")

  # Compute NA overview if not done yet
  if (is.null(na_overview)) {
    na_overview <- .get_NA_overview(se, "log2")
  }

  # =========================================================================
  # RETURN
  # =========================================================================

  list(
    se = se,
    filter_summary = filtered$summary,
    na_overview = na_overview
  )
}
