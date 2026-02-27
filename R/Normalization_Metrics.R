# =============================================================================
# Normalization Quality Metrics
# =============================================================================
#
# Functions for evaluating and comparing proteomics normalization methods:
#   - import_norm_matrices() : Load normalized TSV files into SummarizedExperiment
#   - normalization_metrics(): Orchestrator returning a list of ggplot2 plots
#
# Individual plot functions (10):
#   nm_plot_boxplot, nm_plot_density, nm_plot_pcv,
#   nm_plot_pmad, nm_plot_pev, nm_plot_pca, nm_plot_correlation,
#   nm_plot_mds, nm_plot_scatter, nm_plot_qq
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

#' Per-protein CV averaged across groups (PRONE-style PCV)
#'
#' For each protein, computes CV = 100 * SD / |mean| within each group,
#' then averages across groups. Returns one value per protein.
#'
#' @param mat Numeric matrix (proteins x samples), log2 scale
#' @param groups Factor or character vector of group labels (length = ncol(mat))
#' @return Named numeric vector, one value per protein
#' @keywords internal
.nm_pcv <- function(mat, groups) {
  groups <- as.factor(groups)
  cv_mat <- sapply(levels(groups), function(g) {
    sub <- mat[, groups == g, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      m <- mean(x)
      if (abs(m) < 1e-10) return(NA_real_)
      100 * sd(x) / abs(m)
    })
  })
  if (is.null(dim(cv_mat))) cv_mat else colMeans(cv_mat, na.rm = TRUE)
}

#' Per-protein MAD averaged across groups (PRONE-style PMAD)
#'
#' For each protein, computes MAD = median(|x - median(x)|) within each group,
#' then averages across groups. Returns one value per protein.
#'
#' @param mat Numeric matrix (proteins x samples)
#' @param groups Factor or character vector of group labels
#' @return Named numeric vector, one value per protein
#' @keywords internal
.nm_pmad <- function(mat, groups) {
  groups <- as.factor(groups)
  mad_mat <- sapply(levels(groups), function(g) {
    sub <- mat[, groups == g, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      median(abs(x - median(x)))
    })
  })
  if (is.null(dim(mad_mat))) mad_mat else colMeans(mad_mat, na.rm = TRUE)
}

#' Per-protein variance averaged across groups (PRONE-style PEV)
#'
#' For each protein, computes variance within each group,
#' then averages across groups. Returns one value per protein.
#'
#' @param mat Numeric matrix (proteins x samples)
#' @param groups Factor or character vector of group labels
#' @return Named numeric vector, one value per protein
#' @keywords internal
.nm_pev <- function(mat, groups) {
  groups <- as.factor(groups)
  var_mat <- sapply(levels(groups), function(g) {
    sub <- mat[, groups == g, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      var(x)
    })
  })
  if (is.null(dim(var_mat))) var_mat else colMeans(var_mat, na.rm = TRUE)
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

# Helper: PRONE-style qualitative color palette (RColorBrewer, reversed)
.nm_prone_colors <- function(n) {
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    pal_info   <- RColorBrewer::brewer.pal.info
    qual_pals  <- pal_info[pal_info$category == "qual", ]
    col_vector <- rev(unlist(mapply(RColorBrewer::brewer.pal,
                                    qual_pals$maxcolors,
                                    rownames(qual_pals))))
    rep_len(col_vector, n)
  } else {
    # fallback: evenly spaced hues
    grDevices::hcl(seq(15, 375, length.out = n + 1)[seq_len(n)],
                   l = 65, c = 100)
  }
}

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

# Calcula los límites Tukey de los bigotes por grupo y devuelve un rango
# con padding para usar en coord_cartesian(ylim = ...).
.nm_axis_limits <- function(x, groups, mult = 1.5, padding = 0.05) {
  bounds <- tapply(x, groups, function(vals) {
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) return(c(NA_real_, NA_real_))
    q1  <- stats::quantile(vals, 0.25)
    q3  <- stats::quantile(vals, 0.75)
    iqr <- q3 - q1
    c(max(min(vals), q1 - mult * iqr),
      min(max(vals), q3 + mult * iqr))
  })
  lo  <- min(vapply(bounds, `[`, numeric(1), 1), na.rm = TRUE)
  hi  <- max(vapply(bounds, `[`, numeric(1), 2), na.rm = TRUE)
  pad <- (hi - lo) * padding
  c(lo - pad, hi + pad)
}

# --------------------------------------------------------------------------
# 3. PCV — percentage coefficient of variation
# --------------------------------------------------------------------------

#' PCV boxplot per method (PRONE-style)
#'
#' One boxplot per normalization method showing the per-protein CV distribution
#' (CV averaged across groups for each protein). Lower and less spread values
#' indicate better within-group consistency.
#' With `diff = TRUE`, shows % reduction vs `baseline` as a bar chart.
#'
#' @inheritParams nm_plot_boxplot
#' @param diff Logical. If TRUE, show % reduction vs `baseline`. Default FALSE.
#' @param baseline Character. Assay name used as reference for diff mode.
#'   Default `"log2"`.
#' @return ggplot object.
#' @export
nm_plot_pcv <- function(se, assay_names = NULL,
                        condition_col = "Condition",
                        diff = FALSE, baseline = "log2", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  col_vector  <- .nm_prone_colors(length(assay_names))

  rows <- lapply(assay_names, function(nm) {
    mat <- SummarizedExperiment::assay(se, nm)
    data.frame(Normalization = nm,
               PCV = .nm_pcv(mat, condition),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df$Normalization <- factor(df$Normalization,
                             levels = sort(unique(df$Normalization)))

  if (diff) {
    if (!baseline %in% assay_names)
      stop("baseline '", baseline, "' not found in assay_names.")
    base_mean  <- mean(df$PCV[df$Normalization == baseline], na.rm = TRUE)
    avg_by_nm  <- tapply(df$PCV, df$Normalization, mean, na.rm = TRUE)
    pct_diff   <- (base_mean - avg_by_nm) / abs(base_mean) * 100
    diff_df    <- data.frame(Normalization = names(pct_diff),
                             PCV = as.numeric(pct_diff),
                             stringsAsFactors = FALSE)
    diff_df$Normalization <- factor(diff_df$Normalization,
                                    levels = sort(unique(diff_df$Normalization)))
    ggplot2::ggplot(diff_df,
      ggplot2::aes(x = Normalization, y = PCV, fill = Normalization)) +
      ggplot2::geom_col() +
      ggplot2::geom_label(ggplot2::aes(label = round(PCV, 0)),
                          fill = "white", show.legend = FALSE, color = "black") +
      ggplot2::scale_fill_manual(name = "Normalization Method",
                                 values = col_vector) +
      ggplot2::labs(title = "PCV — % reduction vs baseline",
                    x = "Normalization Method", y = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                          vjust = 0.5))
  } else {
    ggplot2::ggplot(df,
      ggplot2::aes(x = Normalization, y = PCV, fill = Normalization)) +
      ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
      ggplot2::stat_boxplot(geom = "errorbar", width = 0.4, na.rm = TRUE) +
      ggplot2::coord_cartesian(ylim = .nm_axis_limits(df$PCV, df$Normalization)) +
      ggplot2::scale_fill_manual(name = "Normalization Method",
                                 values = col_vector) +
      ggplot2::labs(title = "Intragroup Pooled Coefficient of Variation (PCV)",
                    x = "Normalization Method", y = "PCV") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                          vjust = 0.5))
  }
}

# --------------------------------------------------------------------------
# 4. PMAD — percentage median absolute deviation
# --------------------------------------------------------------------------

#' PMAD boxplot per method (PRONE-style)
#'
#' One boxplot per normalization method showing the per-protein MAD distribution
#' (MAD averaged across groups for each protein). Lower values indicate better
#' within-group consistency.
#' With `diff = TRUE`, shows % reduction vs `baseline` as a bar chart.
#'
#' @inheritParams nm_plot_pcv
#' @return ggplot object.
#' @export
nm_plot_pmad <- function(se, assay_names = NULL,
                         condition_col = "Condition",
                         diff = FALSE, baseline = "log2", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  col_vector  <- .nm_prone_colors(length(assay_names))

  rows <- lapply(assay_names, function(nm) {
    mat <- SummarizedExperiment::assay(se, nm)
    data.frame(Normalization = nm,
               PMAD = .nm_pmad(mat, condition),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df$Normalization <- factor(df$Normalization,
                             levels = sort(unique(df$Normalization)))

  if (diff) {
    if (!baseline %in% assay_names)
      stop("baseline '", baseline, "' not found in assay_names.")
    base_mean <- mean(df$PMAD[df$Normalization == baseline], na.rm = TRUE)
    avg_by_nm <- tapply(df$PMAD, df$Normalization, mean, na.rm = TRUE)
    pct_diff  <- (base_mean - avg_by_nm) / abs(base_mean) * 100
    diff_df   <- data.frame(Normalization = names(pct_diff),
                             PMAD = as.numeric(pct_diff),
                             stringsAsFactors = FALSE)
    diff_df$Normalization <- factor(diff_df$Normalization,
                                    levels = sort(unique(diff_df$Normalization)))
    ggplot2::ggplot(diff_df,
      ggplot2::aes(x = Normalization, y = PMAD, fill = Normalization)) +
      ggplot2::geom_col() +
      ggplot2::geom_label(ggplot2::aes(label = round(PMAD, 0)),
                          fill = "white", show.legend = FALSE, color = "black") +
      ggplot2::scale_fill_manual(name = "Normalization Method",
                                 values = col_vector) +
      ggplot2::labs(title = "PMAD — % reduction vs baseline",
                    x = "Normalization Method", y = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                          vjust = 0.5))
  } else {
    ggplot2::ggplot(df,
      ggplot2::aes(x = Normalization, y = PMAD, fill = Normalization)) +
      ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
      ggplot2::stat_boxplot(geom = "errorbar", width = 0.4, na.rm = TRUE) +
      ggplot2::coord_cartesian(ylim = .nm_axis_limits(df$PMAD, df$Normalization)) +
      ggplot2::scale_fill_manual(name = "Normalization Method",
                                 values = col_vector) +
      ggplot2::labs(title = "Intragroup Pooled Median Absolute Deviation (PMAD)",
                    x = "Normalization Method", y = "PMAD") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                          vjust = 0.5))
  }
}

# --------------------------------------------------------------------------
# 5. PEV — percentage explained variance
# --------------------------------------------------------------------------

#' PEV boxplot per method (PRONE-style)
#'
#' One boxplot per normalization method showing the per-protein variance
#' distribution (variance averaged across groups for each protein).
#' Lower values indicate better within-group consistency.
#' With `diff = TRUE`, shows % reduction vs `baseline` as a bar chart.
#'
#' @inheritParams nm_plot_pcv
#' @return ggplot object.
#' @export
nm_plot_pev <- function(se, assay_names = NULL,
                        condition_col = "Condition",
                        diff = FALSE, baseline = "log2", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  col_vector  <- .nm_prone_colors(length(assay_names))

  rows <- lapply(assay_names, function(nm) {
    mat <- SummarizedExperiment::assay(se, nm)
    data.frame(Normalization = nm,
               PEV = .nm_pev(mat, condition),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df$Normalization <- factor(df$Normalization,
                             levels = sort(unique(df$Normalization)))

  if (diff) {
    if (!baseline %in% assay_names)
      stop("baseline '", baseline, "' not found in assay_names.")
    base_mean <- mean(df$PEV[df$Normalization == baseline], na.rm = TRUE)
    avg_by_nm <- tapply(df$PEV, df$Normalization, mean, na.rm = TRUE)
    pct_diff  <- (base_mean - avg_by_nm) / abs(base_mean) * 100
    diff_df   <- data.frame(Normalization = names(pct_diff),
                             PEV = as.numeric(pct_diff),
                             stringsAsFactors = FALSE)
    diff_df$Normalization <- factor(diff_df$Normalization,
                                    levels = sort(unique(diff_df$Normalization)))
    ggplot2::ggplot(diff_df,
      ggplot2::aes(x = Normalization, y = PEV, fill = Normalization)) +
      ggplot2::geom_col() +
      ggplot2::geom_label(ggplot2::aes(label = round(PEV, 0)),
                          fill = "white", show.legend = FALSE, color = "black") +
      ggplot2::scale_fill_manual(name = "Normalization Method",
                                 values = col_vector) +
      ggplot2::labs(title = "PEV — % reduction vs baseline",
                    x = "Normalization Method", y = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                          vjust = 0.5))
  } else {
    ggplot2::ggplot(df,
      ggplot2::aes(x = Normalization, y = PEV, fill = Normalization)) +
      ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
      ggplot2::stat_boxplot(geom = "errorbar", width = 0.4, na.rm = TRUE) +
      ggplot2::coord_cartesian(ylim = .nm_axis_limits(df$PEV, df$Normalization)) +
      ggplot2::scale_fill_manual(name = "Normalization Method",
                                 values = col_vector) +
      ggplot2::labs(title = "Intragroup Pooled Estimate of Variance (PEV)",
                    x = "Normalization Method", y = "PEV") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                          vjust = 0.5))
  }
}

# --------------------------------------------------------------------------
# 6. PCA — principal component analysis
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
# 7. Correlation — intra-group correlation distribution
# --------------------------------------------------------------------------

#' Intra-group correlation boxplot per method (PRONE-style)
#'
#' Pairwise within-group correlations for each method shown as boxplots
#' with error bars. Higher and less variable values indicate better normalization.
#'
#' @inheritParams nm_plot_boxplot
#' @param cor_method Correlation method passed to `cor()`:
#'   `"pearson"`, `"spearman"`, or `"kendall"`. Default `"pearson"`.
#' @return ggplot object.
#' @export
nm_plot_correlation <- function(se, assay_names = NULL,
                                condition_col = "Condition",
                                cor_method = "pearson", ...) {
  stopifnot(cor_method %in% c("pearson", "spearman", "kendall"))
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)
  col_vector  <- .nm_prone_colors(length(assay_names))

  rows <- lapply(assay_names, function(nm) {
    mat  <- SummarizedExperiment::assay(se, nm)
    cors <- .nm_intragroup_cor(mat, condition, method = cor_method)
    if (length(cors) == 0) return(NULL)
    data.frame(Normalization = nm, Correlation = cors,
               stringsAsFactors = FALSE)
  })
  cor_df <- do.call(rbind, Filter(Negate(is.null), rows))
  cor_df$Normalization <- factor(cor_df$Normalization,
                                 levels = sort(unique(cor_df$Normalization)))

  ggplot2::ggplot(cor_df,
    ggplot2::aes(x = Normalization, y = Correlation, fill = Normalization)) +
    ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.4, na.rm = TRUE) +
    ggplot2::coord_cartesian(ylim = .nm_axis_limits(cor_df$Correlation, cor_df$Normalization)) +
    ggplot2::scale_fill_manual(name = "Normalization Method",
                               values = col_vector) +
    ggplot2::labs(title = "Intragroup Pearson Correlation",
                  x = "Normalization Method",
                  y = paste0(tools::toTitleCase(cor_method), " correlation")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                        vjust = 0.5))
}

# --------------------------------------------------------------------------
# 8. MDS — multidimensional scaling
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
# 9. Scatter — sample-vs-sample (NormalyzerDE style)
# --------------------------------------------------------------------------

#' Sample-vs-sample scatter plot per method (NormalyzerDE style)
#'
#' Plots log2-intensity of one sample against another, one point per protein,
#' for each normalization method. A linear fit and the adjusted R² are overlaid.
#' Default comparison uses the first two columns of the SE (as NormalyzerDE does).
#'
#' @inheritParams nm_plot_boxplot
#' @param sample1 Character. Name of the first sample (x-axis).
#'   Default = first column of `se`.
#' @param sample2 Character. Name of the second sample (y-axis).
#'   Default = second column of `se`.
#' @return ggplot object.
#' @export
nm_plot_scatter <- function(se, assay_names = NULL,
                            condition_col = "Condition",
                            sample1 = NULL, sample2 = NULL, ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  samples     <- colnames(se)

  # Default: first two columns
  sample1 <- sample1 %||% samples[1]
  sample2 <- sample2 %||% samples[2]

  if (!sample1 %in% samples)
    stop("sample1 '", sample1, "' not found in colnames(se).")
  if (!sample2 %in% samples)
    stop("sample2 '", sample2, "' not found in colnames(se).")

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat <- SummarizedExperiment::assay(se, assay_names[i])
    x   <- mat[, sample1]
    y   <- mat[, sample2]
    ok  <- is.finite(x) & is.finite(y)
    r2  <- if (sum(ok) >= 2) {
      summary(lm(y[ok] ~ x[ok]))$adj.r.squared
    } else {
      NA_real_
    }
    df <- data.frame(
      x      = x,
      y      = y,
      Method = assay_names[i],
      R2     = r2,
      stringsAsFactors = FALSE
    )
    rows[[i]] <- df
  }
  scat_df <- do.call(rbind, rows)
  scat_df <- .nm_method_factor(scat_df, assay_names)

  # R² annotation per facet
  r2_df <- unique(scat_df[, c("Method", "R2")])
  r2_df$label <- ifelse(is.na(r2_df$R2), "",
                        paste0("R\u00b2 = ", round(r2_df$R2, 3)))
  r2_df$Method <- factor(r2_df$Method, levels = assay_names)

  # Global axis range for equal scales across facets
  val_range <- range(c(scat_df$x, scat_df$y), na.rm = TRUE)
  x_pos     <- val_range[1] + 0.02 * diff(val_range)
  y_pos     <- val_range[2] - 0.04 * diff(val_range)

  ggplot2::ggplot(scat_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.4, color = "steelblue",
                        na.rm = TRUE) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "firebrick",
                         linewidth = 0.7, na.rm = TRUE, formula = y ~ x) +
    ggplot2::geom_text(data = r2_df,
                       ggplot2::aes(x = x_pos, y = y_pos, label = label),
                       hjust = 0, vjust = 1, size = 3, color = "black",
                       inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ Method, ncol = 2) +
    ggplot2::labs(title = paste0("Sample scatter: ", sample1, " vs ", sample2),
                  x = paste0("log2 Intensity (", sample1, ")"),
                  y = paste0("log2 Intensity (", sample2, ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}

# --------------------------------------------------------------------------
# 10. Q-Q — quantile-quantile normality check (NormalyzerDE style)
# --------------------------------------------------------------------------

#' Q-Q plot per method (NormalyzerDE style)
#'
#' Quantile-quantile plot against a normal distribution for a single sample,
#' one facet per normalization method. Default uses the first column (as
#' NormalyzerDE does).
#'
#' @inheritParams nm_plot_boxplot
#' @param which_sample Character. Name of the sample to inspect.
#'   Default = first column of `se`.
#' @return ggplot object.
#' @export
nm_plot_qq <- function(se, assay_names = NULL,
                       condition_col = "Condition",
                       which_sample = NULL, ...) {
  assay_names  <- .nm_assay_names(se, assay_names)
  samples      <- colnames(se)
  which_sample <- which_sample %||% samples[1]

  if (!which_sample %in% samples)
    stop("which_sample '", which_sample, "' not found in colnames(se).")

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat <- SummarizedExperiment::assay(se, assay_names[i])
    vals <- mat[, which_sample]
    df <- data.frame(
      Value  = vals,
      Method = assay_names[i],
      stringsAsFactors = FALSE
    )
    rows[[i]] <- df
  }
  qq_df <- do.call(rbind, rows)
  qq_df <- .nm_method_factor(qq_df, assay_names)

  ggplot2::ggplot(qq_df, ggplot2::aes(sample = Value)) +
    ggplot2::stat_qq(na.rm = TRUE, size = 0.5, alpha = 0.5,
                     color = "steelblue") +
    ggplot2::stat_qq_line(na.rm = TRUE, color = "firebrick",
                          linewidth = 0.7) +
    ggplot2::facet_wrap(~ Method, ncol = 2) +
    ggplot2::labs(title = paste0("Q-Q plot: sample '", which_sample, "'"),
                  x = "Theoretical quantiles",
                  y = "Sample quantiles") +
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
#'   Valid names: `"boxplot"`, `"density"`, `"pcv"`, `"pmad"`, `"pev"`,
#'   `"pca"`, `"correlation"`, `"mds"`, `"scatter"`, `"qq"`.
#' @param cor_method Correlation method for `nm_plot_correlation()`.
#'   Default `"pearson"`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @return Named list of ggplot objects (or NULL for failed plots).
#'
#' @examples
#' \dontrun{
#' se_nm  <- import_norm_matrices("./results", "./data/metadata.tsv")
#' plots  <- normalization_metrics(se_nm)
#' plots$scatter
#' plots$qq
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
  all_plot_names <- c("boxplot", "density", "pcv", "pmad", "pev",
                      "pca", "correlation", "mds", "scatter", "qq")

  plot_fns <- list(
    boxplot     = function() nm_plot_boxplot(se, assay_names, condition_col),
    density     = function() nm_plot_density(se, assay_names, condition_col),
    pcv         = function() nm_plot_pcv(se, assay_names, condition_col),
    pmad        = function() nm_plot_pmad(se, assay_names, condition_col),
    pev         = function() nm_plot_pev(se, assay_names, condition_col),
    pca         = function() nm_plot_pca(se, assay_names, condition_col),
    correlation = function() nm_plot_correlation(se, assay_names, condition_col,
                                                 cor_method = cor_method),
    mds         = function() nm_plot_mds(se, assay_names, condition_col),
    scatter     = function() nm_plot_scatter(se, assay_names, condition_col),
    qq          = function() nm_plot_qq(se, assay_names, condition_col)
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


  # ---- 2. Generate all 10 quality plots at once ------------------------------

  plots <- normalization_metrics(se_nm)

  # Names of available plots
  names(plots)  # boxplot density pcv pmad pev pca correlation mds scatter qq


  # ---- 3. Inspect individual plots -------------------------------------------

  plots$boxplot     # intensity distribution per sample
  plots$density     # KDE curves per sample
  plots$pcv         # mean CV per condition and method
  plots$pmad        # mean MAD per condition and method
  plots$pev         # mean variance per condition and method
  plots$pca         # PC1 vs PC2, colored by condition
  plots$correlation # intra-group Pearson correlation violin
  plots$mds         # MDS 2D scatter
  plots$scatter     # sample-vs-sample scatter with R²
  plots$qq          # Q-Q normality plot for first sample


  # ---- 4. Single assay, single plot ------------------------------------------

  nm_plot_density(se_nm, assay_names = "cycloess")

  # Scatter between specific samples
  nm_plot_scatter(se_nm, sample1 = "A_1", sample2 = "A_2")

  # Q-Q for a specific sample
  nm_plot_qq(se_nm, which_sample = "B_1")


  # ---- 5. Selective execution via orchestrator --------------------------------

  # Only scatter, Q-Q and PCA for two methods
  subset_plots <- normalization_metrics(
    se_nm,
    assay_names = c("cycloess", "Quantile"),
    plots       = c("scatter", "qq", "pca")
  )
  subset_plots$scatter


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
