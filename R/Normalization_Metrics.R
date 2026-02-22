# =============================================================================
# Normalization Quality Metrics
# =============================================================================
#
# Functions for evaluating and comparing proteomics normalization methods:
#   - import_norm_matrices() : Load normalized TSV files into SummarizedExperiment
#   - normalization_metrics(): Orchestrator returning a list of ggplot2 plots
#
# Individual plot functions (13):
#   nm_plot_boxplot, nm_plot_density, nm_plot_rle, nm_plot_pcv,
#   nm_plot_pmad, nm_plot_pev, nm_plot_pca, nm_plot_correlation,
#   nm_plot_mds, nm_plot_dendrogram, nm_plot_ma, nm_plot_meansd,
#   nm_plot_cv_intensity
#
# References: proteoDA, PRONE, NormalizerDE
#
# Dependencies (required):
#   - ggplot2, dplyr, tidyr
#   - SummarizedExperiment, S4Vectors
#
# Dependencies (optional):
#   - readr    : fast TSV reading (fallback: read.delim)
#   - ggdendro : dendrogram via ggplot2 (fallback: base R dendrogram)
#   - RColorBrewer / viridis : palettes (fallback: ggplot2 defaults)
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# --- Null coalescing operator ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# =============================================================================
# SECTION 1: INTERNAL METRIC HELPERS
# =============================================================================

#' Reshape SE assays to long data.frame for ggplot2
#'
#' @param se SummarizedExperiment object
#' @param assay_names Character vector of assay names to include. NULL = all.
#' @param condition_col Column name in colData with condition labels.
#' @return data.frame with columns: Method, Sample, Condition, Protein, Value
#' @keywords internal
.nm_to_long <- function(se, assay_names = NULL, condition_col = "Condition") {
  if (is.null(assay_names)) assay_names <- SummarizedExperiment::assayNames(se)

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  condition <- if (condition_col %in% colnames(cd)) cd[[condition_col]] else rep("unknown", ncol(se))
  samples   <- colnames(se)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat <- SummarizedExperiment::assay(se, assay_names[i])
    df  <- as.data.frame(mat, stringsAsFactors = FALSE)
    df$Protein <- rownames(mat)
    df_long <- tidyr::pivot_longer(df,
                                   cols      = -Protein,
                                   names_to  = "Sample",
                                   values_to = "Value")
    df_long$Method    <- assay_names[i]
    df_long$Condition <- condition[match(df_long$Sample, samples)]
    rows[[i]] <- df_long
  }

  do.call(rbind, rows)
}

#' Percentage coefficient of variation (PCV) per group
#'
#' For each protein, computes CV = 100 * SD / mean within each group.
#' Returns mean CV across proteins per group.
#'
#' @param mat Numeric matrix (proteins x samples), log2 scale
#' @param groups Factor or character vector of group labels (length = ncol(mat))
#' @return data.frame with columns: Group, PCV
#' @keywords internal
.nm_pcv <- function(mat, groups) {
  groups <- as.factor(groups)
  result <- lapply(levels(groups), function(g) {
    sub_mat <- mat[, groups == g, drop = FALSE]
    # CV per protein: 100 * SD / |mean|
    cv_vals <- apply(sub_mat, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      m <- mean(x)
      if (abs(m) < 1e-10) return(NA_real_)
      100 * sd(x) / abs(m)
    })
    data.frame(Group = g, PCV = mean(cv_vals, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}

#' Percentage median absolute deviation (PMAD) per group
#'
#' For each protein, computes MAD = median(|x - median(x)|) within each group.
#' Returns mean MAD across proteins per group.
#'
#' @param mat Numeric matrix (proteins x samples)
#' @param groups Factor or character vector of group labels
#' @return data.frame with columns: Group, PMAD
#' @keywords internal
.nm_pmad <- function(mat, groups) {
  groups <- as.factor(groups)
  result <- lapply(levels(groups), function(g) {
    sub_mat <- mat[, groups == g, drop = FALSE]
    mad_vals <- apply(sub_mat, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      median(abs(x - median(x)))
    })
    data.frame(Group = g, PMAD = mean(mad_vals, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}

#' Percentage explained variance (PEV) per group
#'
#' For each protein, computes variance within each group.
#' Returns mean variance across proteins per group.
#'
#' @param mat Numeric matrix (proteins x samples)
#' @param groups Factor or character vector of group labels
#' @return data.frame with columns: Group, PEV
#' @keywords internal
.nm_pev <- function(mat, groups) {
  groups <- as.factor(groups)
  result <- lapply(levels(groups), function(g) {
    sub_mat <- mat[, groups == g, drop = FALSE]
    var_vals <- apply(sub_mat, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      var(x)
    })
    data.frame(Group = g, PEV = mean(var_vals, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}

#' Relative log expression (RLE) matrix
#'
#' Computes row-wise deviation from each protein's median.
#'
#' @param mat Numeric matrix (proteins x samples)
#' @return Numeric matrix of same dimensions: value - row median
#' @keywords internal
.nm_rle <- function(mat) {
  row_med <- apply(mat, 1, median, na.rm = TRUE)
  mat - row_med
}

#' Intra-group Pearson correlations
#'
#' Computes pairwise Pearson correlations between all sample pairs
#' within each group, using complete observations only.
#'
#' @param mat Numeric matrix (proteins x samples)
#' @param groups Factor or character vector of group labels
#' @return Numeric vector of correlation values
#' @keywords internal
.nm_intragroup_cor <- function(mat, groups, method = "pearson") {
  groups  <- as.factor(groups)
  cor_all <- c()
  for (g in levels(groups)) {
    idx <- which(groups == g)
    if (length(idx) < 2) next
    sub  <- mat[, idx, drop = FALSE]
    keep <- complete.cases(sub)
    if (sum(keep) < 2) next
    cor_mat <- cor(sub[keep, ], method = method, use = "pairwise.complete.obs")
    upper   <- cor_mat[upper.tri(cor_mat)]
    cor_all <- c(cor_all, upper)
  }
  cor_all
}

# =============================================================================
# SECTION 2: IMPORT FUNCTION
# =============================================================================

#' Import normalized matrices into a SummarizedExperiment
#'
#' Reads all TSV files matching a pattern from a directory, one assay per file.
#' Each TSV must have column 1 named `ProteinGroups` and remaining columns as
#' sample names matching the `sample_col` in the metadata file.
#'
#' @param tsv_dir Character. Directory containing the normalized TSV files.
#' @param metadata_path Character. Path to a TSV or CSV file with at minimum
#'   columns `sample_col` and `condition_col`.
#' @param pattern Character. Regex to filter files in `tsv_dir`.
#'   Default: `"matrix_log2_.*\\.tsv$"`.
#' @param method_names Character vector. Assay names to assign (in order).
#'   NULL (default) = extracted from file names between `matrix_log2_` and `.tsv`.
#' @param condition_col Character. Column in metadata with condition/group labels.
#' @param sample_col Character. Column in metadata with sample names.
#' @return A SummarizedExperiment with one assay per TSV, colData from metadata,
#'   and rowData containing `Protein.IDs`.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices(
#'   tsv_dir       = "./results",
#'   metadata_path = "./data/metadata.tsv"
#' )
#' SummarizedExperiment::assayNames(se_nm)
#' }
#' @export
import_norm_matrices <- function(tsv_dir,
                                 metadata_path,
                                 pattern       = "matrix_log2_.*\\.tsv$",
                                 method_names  = NULL,
                                 condition_col = "Condition",
                                 sample_col    = "Column") {
  # --- Required packages ---
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("Package 'SummarizedExperiment' is required.")
  if (!requireNamespace("S4Vectors", quietly = TRUE))
    stop("Package 'S4Vectors' is required.")

  # --- Find TSV files ---
  files <- list.files(tsv_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0)
    stop("No files matching pattern '", pattern, "' found in: ", tsv_dir)

  # --- Determine assay names ---
  if (is.null(method_names)) {
    base_names   <- basename(files)
    method_names <- sub("^matrix_log2_", "", base_names)
    method_names <- sub("\\.tsv$", "", method_names)
  }
  if (length(method_names) != length(files))
    stop("Length of method_names (", length(method_names),
         ") must match number of files (", length(files), ").")

  # --- Read metadata ---
  meta_ext <- tolower(tools::file_ext(metadata_path))
  if (requireNamespace("readr", quietly = TRUE)) {
    meta <- if (meta_ext == "csv") {
      readr::read_csv(metadata_path, show_col_types = FALSE)
    } else {
      readr::read_tsv(metadata_path, show_col_types = FALSE)
    }
    meta <- as.data.frame(meta)
  } else {
    sep  <- if (meta_ext == "csv") "," else "\t"
    meta <- read.delim(metadata_path, sep = sep, stringsAsFactors = FALSE,
                       check.names = FALSE)
  }

  if (!sample_col %in% colnames(meta))
    stop("Column '", sample_col, "' not found in metadata.")
  if (!condition_col %in% colnames(meta))
    stop("Column '", condition_col, "' not found in metadata.")

  rownames(meta) <- meta[[sample_col]]

  # --- Read first file to get sample names and protein IDs ---
  .read_tsv_file <- function(path) {
    if (requireNamespace("readr", quietly = TRUE)) {
      as.data.frame(readr::read_tsv(path, show_col_types = FALSE),
                    stringsAsFactors = FALSE)
    } else {
      read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    }
  }

  first_df  <- .read_tsv_file(files[1])
  protein_col <- colnames(first_df)[1]  # should be "ProteinGroups"
  proteins  <- first_df[[protein_col]]
  sample_cols <- setdiff(colnames(first_df), protein_col)

  # --- Validate samples against metadata ---
  missing_meta <- setdiff(sample_cols, meta[[sample_col]])
  if (length(missing_meta) > 0)
    stop("Samples in TSV not found in metadata: ",
         paste(missing_meta, collapse = ", "))

  # Align metadata to sample order in TSV
  meta_aligned       <- meta[sample_cols, , drop = FALSE]
  rownames(meta_aligned) <- sample_cols

  # --- Build assay list ---
  assay_list <- vector("list", length(files))
  names(assay_list) <- method_names

  for (i in seq_along(files)) {
    df    <- .read_tsv_file(files[i])
    mat   <- as.matrix(df[, sample_cols, drop = FALSE])
    mode(mat) <- "numeric"
    rownames(mat) <- proteins
    assay_list[[i]] <- mat
  }

  # --- Build SummarizedExperiment ---
  row_data <- S4Vectors::DataFrame(Protein.IDs = proteins)
  rownames(row_data) <- proteins

  col_data <- S4Vectors::DataFrame(meta_aligned)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = assay_list,
    colData = col_data,
    rowData = row_data
  )

  message("import_norm_matrices: loaded ", length(files), " assay(s) — ",
          paste(method_names, collapse = ", "))
  message("  Proteins: ", nrow(se), " | Samples: ", ncol(se))
  se
}

# =============================================================================
# SECTION 3: INDIVIDUAL PLOT FUNCTIONS
# =============================================================================

# Helper: resolve assay_names (NULL → all)
.nm_assay_names <- function(se, assay_names) {
  assay_names %||% SummarizedExperiment::assayNames(se)
}

# Helper: get condition vector aligned to colnames(se)
.nm_condition <- function(se, condition_col) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if (condition_col %in% colnames(cd)) cd[[condition_col]] else rep("unknown", ncol(se))
}

# Helper: order Method factor by assay order
.nm_method_factor <- function(df, assay_names) {
  df$Method <- factor(df$Method, levels = assay_names)
  df
}

# --------------------------------------------------------------------------
# 1. Boxplot — intensity distribution per sample
# --------------------------------------------------------------------------

#' Intensity boxplot per method
#'
#' Displays distribution of protein intensities for each sample,
#' faceted by normalization method.
#'
#' @param se SummarizedExperiment from `import_norm_matrices()`.
#' @param assay_names Character vector of assay names to include. NULL = all.
#' @param condition_col Column in colData with condition labels.
#' @param ... Additional arguments (unused).
#' @return ggplot object.
#' @export
nm_plot_boxplot <- function(se, assay_names = NULL,
                             condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  long_df     <- .nm_to_long(se, assay_names, condition_col)
  long_df     <- .nm_method_factor(long_df, assay_names)

  ggplot2::ggplot(long_df,
    ggplot2::aes(x = Sample, y = Value, fill = Condition)) +
    ggplot2::geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4,
                          na.rm = TRUE) +
    ggplot2::facet_wrap(~ Method, ncol = 2, scales = "free_x") +
    ggplot2::labs(title = "Intensity distribution by sample",
                  x = NULL, y = "log2 Intensity") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                        size = 7),
                   strip.text  = ggplot2::element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# 2. Density — KDE curves per sample
# --------------------------------------------------------------------------

#' Density plot per method
#'
#' KDE curves for each sample's intensity distribution, faceted by method.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_density <- function(se, assay_names = NULL,
                             condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  long_df     <- .nm_to_long(se, assay_names, condition_col)
  long_df     <- .nm_method_factor(long_df, assay_names)

  ggplot2::ggplot(long_df,
    ggplot2::aes(x = Value, group = Sample, color = Condition)) +
    ggplot2::geom_density(na.rm = TRUE, linewidth = 0.5, alpha = 0.8) +
    ggplot2::facet_wrap(~ Method, ncol = 2) +
    ggplot2::labs(title = "Intensity density per sample",
                  x = "log2 Intensity", y = "Density") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# 3. RLE — relative log expression
# --------------------------------------------------------------------------

#' RLE boxplot per method
#'
#' Relative log expression (value - row median) per sample,
#' faceted by method. Boxes should be centered at y = 0 for ideal normalization.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_rle <- function(se, assay_names = NULL,
                        condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  samples     <- colnames(se)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat  <- SummarizedExperiment::assay(se, assay_names[i])
    rle  <- .nm_rle(mat)
    df   <- as.data.frame(rle, stringsAsFactors = FALSE)
    df$Protein <- rownames(rle)
    df_long <- tidyr::pivot_longer(df, cols = -Protein,
                                   names_to = "Sample", values_to = "RLE")
    df_long$Method    <- assay_names[i]
    df_long$Condition <- condition[match(df_long$Sample, samples)]
    rows[[i]] <- df_long
  }
  rle_df <- do.call(rbind, rows)
  rle_df <- .nm_method_factor(rle_df, assay_names)

  ggplot2::ggplot(rle_df,
    ggplot2::aes(x = Sample, y = RLE, fill = Condition)) +
    ggplot2::geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3,
                          na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "firebrick", linewidth = 0.6) +
    ggplot2::facet_wrap(~ Method, ncol = 2, scales = "free_x") +
    ggplot2::labs(title = "Relative Log Expression (RLE)",
                  x = NULL, y = "RLE (value - row median)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                        size = 7),
                   strip.text  = ggplot2::element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# 4. PCV — percentage coefficient of variation
# --------------------------------------------------------------------------

#' PCV point-range plot per method
#'
#' Mean coefficient of variation (%) per condition across all methods.
#' Lower values indicate better within-group consistency.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_pcv <- function(se, assay_names = NULL,
                        condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat      <- SummarizedExperiment::assay(se, assay_names[i])
    pcv_df   <- .nm_pcv(mat, condition)
    pcv_df$Method <- assay_names[i]
    rows[[i]] <- pcv_df
  }
  df <- do.call(rbind, rows)
  df <- .nm_method_factor(df, assay_names)

  ggplot2::ggplot(df,
    ggplot2::aes(x = Method, y = PCV, color = Group, group = Group)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(title = "Percentage Coefficient of Variation (PCV)",
                  x = "Normalization method", y = "Mean PCV (%)",
                  color = "Condition") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# --------------------------------------------------------------------------
# 5. PMAD — percentage median absolute deviation
# --------------------------------------------------------------------------

#' PMAD point-range plot per method
#'
#' Mean median absolute deviation (log2) per condition across methods.
#' Lower values indicate better within-group consistency.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_pmad <- function(se, assay_names = NULL,
                         condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat     <- SummarizedExperiment::assay(se, assay_names[i])
    pmad_df <- .nm_pmad(mat, condition)
    pmad_df$Method <- assay_names[i]
    rows[[i]] <- pmad_df
  }
  df <- do.call(rbind, rows)
  df <- .nm_method_factor(df, assay_names)

  ggplot2::ggplot(df,
    ggplot2::aes(x = Method, y = PMAD, color = Group, group = Group)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(title = "Percentage Median Absolute Deviation (PMAD)",
                  x = "Normalization method", y = "Mean PMAD",
                  color = "Condition") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# --------------------------------------------------------------------------
# 6. PEV — percentage explained variance
# --------------------------------------------------------------------------

#' PEV point-range plot per method
#'
#' Mean within-group variance per condition across methods.
#' Lower values indicate better within-group consistency.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_pev <- function(se, assay_names = NULL,
                        condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat    <- SummarizedExperiment::assay(se, assay_names[i])
    pev_df <- .nm_pev(mat, condition)
    pev_df$Method <- assay_names[i]
    rows[[i]] <- pev_df
  }
  df <- do.call(rbind, rows)
  df <- .nm_method_factor(df, assay_names)

  ggplot2::ggplot(df,
    ggplot2::aes(x = Method, y = PEV, color = Group, group = Group)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(title = "Percentage Explained Variance (PEV)",
                  x = "Normalization method", y = "Mean variance",
                  color = "Condition") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# --------------------------------------------------------------------------
# 7. PCA — principal component analysis
# --------------------------------------------------------------------------

#' PCA scatter plot per method
#'
#' PC1 vs PC2 scatter, colored by condition, faceted by method.
#' Percentage of variance explained shown on each axis.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_pca <- function(se, assay_names = NULL,
                        condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  samples     <- colnames(se)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat     <- SummarizedExperiment::assay(se, assay_names[i])
    # Remove rows with any NA then transpose
    mat_ok  <- mat[complete.cases(mat), ]
    if (nrow(mat_ok) < 2) {
      message("nm_plot_pca: insufficient complete rows for '", assay_names[i],
              "', skipping.")
      next
    }
    pca_res <- prcomp(t(mat_ok), scale. = FALSE, center = TRUE)
    pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
    df      <- data.frame(
      PC1       = pca_res$x[, 1],
      PC2       = pca_res$x[, 2],
      Sample    = rownames(pca_res$x),
      Condition = condition[match(rownames(pca_res$x), samples)],
      Method    = assay_names[i],
      Pct1      = pct_var[1],
      Pct2      = pct_var[2],
      stringsAsFactors = FALSE
    )
    rows[[i]] <- df
  }
  pca_df <- do.call(rbind, rows)
  if (is.null(pca_df) || nrow(pca_df) == 0)
    stop("nm_plot_pca: no valid data for PCA.")

  # Build per-method axis labels (use first occurrence of Pct1/Pct2)
  labels_df <- unique(pca_df[, c("Method", "Pct1", "Pct2")])
  pca_df    <- .nm_method_factor(pca_df, assay_names)

  # One shared facet label is not ideal when % var differs; use subtitles via
  # labeller (simple approach: embed % in facet label)
  pca_df$Facet <- paste0(pca_df$Method,
                         " (PC1=", pca_df$Pct1, "%, PC2=", pca_df$Pct2, "%)")

  ggplot2::ggplot(pca_df,
    ggplot2::aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~ Facet, ncol = 2, scales = "free") +
    ggplot2::labs(title = "PCA — PC1 vs PC2", x = "PC1", y = "PC2") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8))
}

# --------------------------------------------------------------------------
# 8. Correlation — intra-group correlation distribution
# --------------------------------------------------------------------------

#' Intra-group correlation violin per method
#'
#' Distribution of pairwise within-group Pearson correlations for each method.
#' Higher and less variable values indicate better normalization.
#'
#' @inheritParams nm_plot_boxplot
#' @param cor_method Correlation method passed to `cor()`. Default "pearson".
#' @return ggplot object.
#' @export
nm_plot_correlation <- function(se, assay_names = NULL,
                                condition_col = "Condition",
                                cor_method = "pearson", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat    <- SummarizedExperiment::assay(se, assay_names[i])
    cors   <- .nm_intragroup_cor(mat, condition, method = cor_method)
    if (length(cors) == 0) next
    rows[[i]] <- data.frame(Method = assay_names[i], Correlation = cors,
                            stringsAsFactors = FALSE)
  }
  cor_df <- do.call(rbind, rows)
  cor_df <- .nm_method_factor(cor_df, assay_names)

  ggplot2::ggplot(cor_df,
    ggplot2::aes(x = Method, y = Correlation, fill = Method)) +
    ggplot2::geom_violin(trim = FALSE, na.rm = TRUE, alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.08, outlier.size = 0.6,
                          fill = "white", na.rm = TRUE) +
    ggplot2::labs(title = paste0("Intra-group ", cor_method,
                                 " correlation distribution"),
                  x = "Normalization method", y = "Pearson r") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "none")
}

# --------------------------------------------------------------------------
# 9. MDS — multidimensional scaling
# --------------------------------------------------------------------------

#' MDS 2D scatter per method
#'
#' Classical MDS on Euclidean distances of scaled, complete-case data.
#' Faceted by method, colored by condition.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_mds <- function(se, assay_names = NULL,
                        condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  samples     <- colnames(se)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat    <- SummarizedExperiment::assay(se, assay_names[i])
    mat_ok <- mat[complete.cases(mat), ]
    if (nrow(mat_ok) < 2) next
    d      <- dist(scale(t(mat_ok)))
    mds    <- cmdscale(d, k = 2)
    df     <- data.frame(
      MDS1      = mds[, 1],
      MDS2      = mds[, 2],
      Sample    = rownames(mds),
      Condition = condition[match(rownames(mds), samples)],
      Method    = assay_names[i],
      stringsAsFactors = FALSE
    )
    rows[[i]] <- df
  }
  mds_df <- do.call(rbind, rows)
  mds_df <- .nm_method_factor(mds_df, assay_names)

  ggplot2::ggplot(mds_df,
    ggplot2::aes(x = MDS1, y = MDS2, color = Condition, label = Sample)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~ Method, ncol = 2, scales = "free") +
    ggplot2::labs(title = "Multidimensional Scaling (MDS)",
                  x = "MDS1", y = "MDS2") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# 10. Dendrogram — hierarchical clustering
# --------------------------------------------------------------------------

#' Hierarchical clustering dendrogram per method
#'
#' Uses `hclust(dist(t(scale(mat))), method = "average")`.
#' If `ggdendro` is available, returns a ggplot; otherwise uses base R graphics.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot (if ggdendro available) or NULL (base R plot rendered).
#' @export
nm_plot_dendrogram <- function(se, assay_names = NULL,
                               condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  samples     <- colnames(se)
  has_ggdendro <- requireNamespace("ggdendro", quietly = TRUE)

  if (!has_ggdendro) {
    # Base R fallback: render all in one panel, warn user
    message("nm_plot_dendrogram: 'ggdendro' not available — ",
            "rendering base R dendrograms (not ggplot2).")
    n_methods <- length(assay_names)
    old_par   <- graphics::par(mfrow = c(ceiling(n_methods / 2), 2))
    on.exit(graphics::par(old_par))
    for (nm in assay_names) {
      mat    <- SummarizedExperiment::assay(se, nm)
      mat_ok <- mat[complete.cases(mat), ]
      if (nrow(mat_ok) < 2) next
      hc  <- stats::hclust(dist(t(scale(mat_ok))), method = "average")
      cond_cols <- condition[match(hc$labels, samples)]
      dend <- stats::as.dendrogram(hc)
      col_map <- setNames(grDevices::rainbow(length(unique(condition))),
                          unique(condition))
      label_cols <- col_map[cond_cols]
      graphics::plot(dend, main = nm, xlab = "", ylab = "Height",
                     nodePar = list(lab.col = label_cols, pch = NA))
    }
    return(invisible(NULL))
  }

  # ggdendro path: build one ggplot with facets
  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat    <- SummarizedExperiment::assay(se, assay_names[i])
    mat_ok <- mat[complete.cases(mat), ]
    if (nrow(mat_ok) < 2) next
    hc      <- stats::hclust(dist(t(scale(mat_ok))), method = "average")
    dend_df <- ggdendro::dendro_data(hc)
    seg_df  <- ggdendro::segment(dend_df)
    lab_df  <- ggdendro::label(dend_df)
    lab_df$Condition <- condition[match(lab_df$label, samples)]
    seg_df$Method   <- assay_names[i]
    lab_df$Method   <- assay_names[i]
    rows[[i]] <- list(seg = seg_df, lab = lab_df)
  }
  rows   <- Filter(Negate(is.null), rows)
  seg_all <- do.call(rbind, lapply(rows, `[[`, "seg"))
  lab_all <- do.call(rbind, lapply(rows, `[[`, "lab"))
  seg_all$Method <- factor(seg_all$Method, levels = assay_names)
  lab_all$Method <- factor(lab_all$Method, levels = assay_names)

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg_all,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::geom_text(data = lab_all,
      ggplot2::aes(x = x, y = -0.01, label = label, color = Condition),
      angle = 90, hjust = 1, size = 2.5) +
    ggplot2::facet_wrap(~ Method, ncol = 2, scales = "free_x") +
    ggplot2::labs(title = "Hierarchical clustering dendrogram",
                  x = NULL, y = "Height") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   strip.text   = ggplot2::element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# 11. MA — M-A plot
# --------------------------------------------------------------------------

#' MA plot per method
#'
#' For each protein in each sample, plots M = (sample - group mean) vs
#' A = group mean intensity. A LOESS smooth and reference y = 0 are added.
#' Faceted by method.
#'
#' @inheritParams nm_plot_boxplot
#' @param max_proteins Integer. Maximum proteins to plot (random sample for speed).
#'   Default 2000.
#' @return ggplot object.
#' @export
nm_plot_ma <- function(se, assay_names = NULL,
                       condition_col = "Condition",
                       max_proteins = 2000L, ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  samples     <- colnames(se)
  groups      <- as.factor(condition)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat    <- SummarizedExperiment::assay(se, assay_names[i])
    # Group means per protein
    group_mean <- sapply(levels(groups), function(g) {
      idx <- which(groups == g)
      rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
    })
    # Long format: one row per (protein, sample)
    ma_rows <- vector("list", ncol(mat))
    for (s in seq_len(ncol(mat))) {
      g     <- as.character(groups[s])
      A_val <- group_mean[, g]
      M_val <- mat[, s] - A_val
      ma_rows[[s]] <- data.frame(A = A_val, M = M_val,
                                 Protein = rownames(mat),
                                 Sample = samples[s],
                                 Method = assay_names[i],
                                 stringsAsFactors = FALSE)
    }
    df <- do.call(rbind, ma_rows)
    df <- df[is.finite(df$A) & is.finite(df$M), ]
    # Subsample for speed
    if (nrow(df) > max_proteins * ncol(mat)) {
      set.seed(42)
      df <- df[sample(nrow(df), min(nrow(df), max_proteins * ncol(mat))), ]
    }
    rows[[i]] <- df
  }
  ma_df <- do.call(rbind, rows)
  ma_df <- .nm_method_factor(ma_df, assay_names)

  ggplot2::ggplot(ma_df, ggplot2::aes(x = A, y = M)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.2, color = "steelblue",
                        na.rm = TRUE) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "firebrick",
                         linewidth = 0.8, na.rm = TRUE, formula = y ~ x) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "black", linewidth = 0.5) +
    ggplot2::facet_wrap(~ Method, ncol = 2) +
    ggplot2::labs(title = "MA plot (M = sample − group mean)",
                  x = "A (group mean intensity)", y = "M (deviation)") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# 12. Mean-SD — SD vs mean
# --------------------------------------------------------------------------

#' Mean-SD plot per method
#'
#' Standard deviation vs mean intensity per protein, with LOESS trend.
#' A flat trend indicates stabilized variance (ideal after normalization).
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_meansd <- function(se, assay_names = NULL,
                           condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat  <- SummarizedExperiment::assay(se, assay_names[i])
    Mean <- rowMeans(mat, na.rm = TRUE)
    SD   <- apply(mat, 1, sd, na.rm = TRUE)
    df   <- data.frame(Mean = Mean, SD = SD, Method = assay_names[i],
                       stringsAsFactors = FALSE)
    df   <- df[is.finite(df$Mean) & is.finite(df$SD), ]
    rows[[i]] <- df
  }
  ms_df <- do.call(rbind, rows)
  ms_df <- .nm_method_factor(ms_df, assay_names)

  ggplot2::ggplot(ms_df, ggplot2::aes(x = Mean, y = SD)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.3, color = "steelblue",
                        na.rm = TRUE) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "firebrick",
                         linewidth = 0.8, na.rm = TRUE, formula = y ~ x) +
    ggplot2::facet_wrap(~ Method, ncol = 2) +
    ggplot2::labs(title = "Mean-SD relationship per protein",
                  x = "Mean (log2 Intensity)", y = "SD") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# 13. CV vs intensity — bonus
# --------------------------------------------------------------------------

#' CV vs mean intensity scatter per method
#'
#' Coefficient of variation (%) vs mean intensity per protein, with LOESS trend.
#' Ideal normalization reduces CV heteroscedasticity.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#' @export
nm_plot_cv_intensity <- function(se, assay_names = NULL,
                                 condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat  <- SummarizedExperiment::assay(se, assay_names[i])
    Mean <- rowMeans(mat, na.rm = TRUE)
    SD   <- apply(mat, 1, sd, na.rm = TRUE)
    CV   <- 100 * SD / abs(Mean)
    df   <- data.frame(Mean = Mean, CV = CV, Method = assay_names[i],
                       stringsAsFactors = FALSE)
    df   <- df[is.finite(df$Mean) & is.finite(df$CV), ]
    rows[[i]] <- df
  }
  cv_df <- do.call(rbind, rows)
  cv_df <- .nm_method_factor(cv_df, assay_names)

  ggplot2::ggplot(cv_df, ggplot2::aes(x = Mean, y = CV)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.3, color = "steelblue",
                        na.rm = TRUE) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "firebrick",
                         linewidth = 0.8, na.rm = TRUE, formula = y ~ x) +
    ggplot2::facet_wrap(~ Method, ncol = 2) +
    ggplot2::labs(title = "CV vs mean intensity per protein",
                  x = "Mean (log2 Intensity)", y = "CV (%)") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}

# =============================================================================
# SECTION 4: MAIN ORCHESTRATOR
# =============================================================================

#' Generate quality metric plots for normalization comparison
#'
#' Calls all (or selected) `nm_plot_*()` functions on a SummarizedExperiment
#' produced by `import_norm_matrices()`. Each plot is wrapped in `tryCatch()`
#' so that a failure in one plot does not abort the entire run.
#'
#' @param se SummarizedExperiment from `import_norm_matrices()`.
#' @param assay_names Character vector of assay names to include.
#'   NULL (default) = all assays.
#' @param condition_col Column name in `colData(se)` with condition labels.
#'   Default `"Condition"`.
#' @param plots Character vector of plot names to generate, or `"all"` (default).
#'   Valid names: `"boxplot"`, `"density"`, `"rle"`, `"pcv"`, `"pmad"`,
#'   `"pev"`, `"pca"`, `"correlation"`, `"mds"`, `"dendrogram"`, `"ma"`,
#'   `"meansd"`, `"cv_intensity"`.
#' @param cor_method Correlation method for `nm_plot_correlation()`.
#'   Default `"pearson"`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @return Named list of ggplot objects (or NULL for failed plots).
#'
#' @examples
#' \dontrun{
#' se_nm  <- import_norm_matrices("./results", "./data/metadata.tsv")
#' plots  <- normalization_metrics(se_nm)
#' plots$rle
#' plots$pca
#' plots$boxplot
#' }
#' @export
normalization_metrics <- function(se,
                                  assay_names   = NULL,
                                  condition_col = "Condition",
                                  plots         = "all",
                                  cor_method    = "pearson",
                                  verbose       = TRUE) {
  # --- Required packages check ---
  for (pkg in c("ggplot2", "dplyr", "tidyr", "SummarizedExperiment", "S4Vectors")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for normalization_metrics().")
  }

  # --- Resolve assay names ---
  all_assays  <- SummarizedExperiment::assayNames(se)
  assay_names <- assay_names %||% all_assays
  missing_a   <- setdiff(assay_names, all_assays)
  if (length(missing_a) > 0)
    stop("Assay(s) not found in SE: ", paste(missing_a, collapse = ", "))

  # --- Plot registry ---
  all_plot_names <- c("boxplot", "density", "rle", "pcv", "pmad", "pev",
                      "pca", "correlation", "mds", "dendrogram", "ma",
                      "meansd", "cv_intensity")

  plot_fns <- list(
    boxplot      = function() nm_plot_boxplot(se, assay_names, condition_col),
    density      = function() nm_plot_density(se, assay_names, condition_col),
    rle          = function() nm_plot_rle(se, assay_names, condition_col),
    pcv          = function() nm_plot_pcv(se, assay_names, condition_col),
    pmad         = function() nm_plot_pmad(se, assay_names, condition_col),
    pev          = function() nm_plot_pev(se, assay_names, condition_col),
    pca          = function() nm_plot_pca(se, assay_names, condition_col),
    correlation  = function() nm_plot_correlation(se, assay_names, condition_col,
                                                  cor_method = cor_method),
    mds          = function() nm_plot_mds(se, assay_names, condition_col),
    dendrogram   = function() nm_plot_dendrogram(se, assay_names, condition_col),
    ma           = function() nm_plot_ma(se, assay_names, condition_col),
    meansd       = function() nm_plot_meansd(se, assay_names, condition_col),
    cv_intensity = function() nm_plot_cv_intensity(se, assay_names, condition_col)
  )

  # --- Determine which plots to run ---
  if (identical(plots, "all")) {
    selected <- all_plot_names
  } else {
    unknown <- setdiff(plots, all_plot_names)
    if (length(unknown) > 0)
      warning("Unknown plot name(s) ignored: ", paste(unknown, collapse = ", "))
    selected <- intersect(plots, all_plot_names)
    if (length(selected) == 0)
      stop("No valid plot names provided.")
  }

  # --- Execute plots with error handling ---
  result <- vector("list", length(selected))
  names(result) <- selected

  for (nm in selected) {
    if (verbose) message("  Generating plot: ", nm, " ...")
    result[[nm]] <- tryCatch(
      plot_fns[[nm]](),
      error = function(e) {
        warning("nm_plot_", nm, "() failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  n_ok   <- sum(!sapply(result, is.null))
  n_fail <- length(selected) - n_ok
  if (verbose) {
    message("normalization_metrics: ", n_ok, " plot(s) generated",
            if (n_fail > 0) paste0(", ", n_fail, " failed") else ".")
  }

  result
}

# =============================================================================
# EXAMPLE WORKFLOW
# =============================================================================
#
# Assumes:
#   - ./results/ contains files like:
#       matrix_log2_cycloess.tsv
#       matrix_log2_Quantile.tsv
#       matrix_log2_vsn.tsv
#   - ./data/metadata.tsv has at least two columns: Column, Condition
#
# Run with:
#   source("R/Normalization_Metrics.R")
# -----------------------------------------------------------------------------

if (FALSE) {

  # ---- 1. Load all normalized matrices into a SummarizedExperiment ----------

  se_nm <- import_norm_matrices(
    tsv_dir       = "./results",
    metadata_path = "./data/metadata.tsv",   # columns: Column, Condition
    pattern       = "matrix_log2_.*\\.tsv$"
  )

  # Check loaded assays and dimensions
  SummarizedExperiment::assayNames(se_nm)  # e.g. "cycloess", "Quantile", "vsn"
  dim(se_nm)                               # proteins x samples


  # ---- 2. Generate all 13 quality plots at once ------------------------------

  plots <- normalization_metrics(se_nm)

  # Names of available plots
  names(plots)


  # ---- 3. Inspect individual plots -------------------------------------------

  plots$boxplot      # intensity distribution per sample
  plots$density      # KDE curves per sample
  plots$rle          # RLE — boxes should be centered at y = 0
  plots$pca          # PC1 vs PC2, colored by condition
  plots$correlation  # intra-group Pearson correlation violin
  plots$mds          # MDS 2D scatter
  plots$dendrogram   # hierarchical clustering
  plots$ma           # MA plot (M = sample − group mean)
  plots$meansd       # SD vs mean — flat trend = ideal
  plots$cv_intensity # CV(%) vs mean intensity
  plots$pcv          # mean CV per condition and method
  plots$pmad         # mean MAD per condition and method
  plots$pev          # mean variance per condition and method


  # ---- 4. Single assay, single plot ------------------------------------------

  nm_plot_density(se_nm, assay_names = "cycloess")

  nm_plot_rle(se_nm, assay_names = c("cycloess", "Quantile"))


  # ---- 5. Selective execution via orchestrator --------------------------------

  # Only RLE, PCA and correlation for two methods
  subset_plots <- normalization_metrics(
    se_nm,
    assay_names = c("cycloess", "Quantile"),
    plots       = c("rle", "pca", "correlation")
  )
  subset_plots$rle


  # ---- 6. Export plots to PNG -------------------------------------------------

  output_dir <- "./results/normalization_metrics"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (plot_name in names(plots)) {
    p <- plots[[plot_name]]
    if (is.null(p)) next
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0("nm_", plot_name, ".png")),
      plot     = p,
      width    = 12,
      height   = 8,
      dpi      = 150
    )
  }

}
