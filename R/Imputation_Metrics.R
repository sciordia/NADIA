# =============================================================================
# Imputation Quality Metrics
# =============================================================================
#
# Functions for evaluating and comparing proteomics imputation methods using
# ground-truth simulation (NAguideR framework):
#   - im_prepare_se()       : Build SE from preprocessing + best normalization
#   - import_imp_matrices() : Load imputed TSV files into SummarizedExperiment
#   - im_compute_metrics()  : Compute NRMSE, SOR, PSS, ACC_OI + ranking
#   - imputation_metrics()  : Orchestrator returning plots + metrics_table
#
# Individual plot functions (6):
#   im_plot_nrmse, im_plot_sor, im_plot_pss,
#   im_plot_acc_oi, im_plot_ranking, im_plot_metrics
#
# References: NAguideR (Wang et al., 2020)
#
# Dependencies (required):
#   - ggplot2, dplyr, tidyr
#   - SummarizedExperiment, S4Vectors
#
# Dependencies (optional):
#   - vegan    : PSS (Procrustes, fallback: NA)
#   - readr    : fast TSV reading (fallback: read.delim)
#   - Imputation.R deps: per-method optional packages
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# --- Null coalescing operator ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# --- Self-dir sourcing for Imputation.R ---
.self_dir <- if (sys.nframe() > 0) dirname(sys.frame(1)$ofile) else "R"

if (!exists(".dispatch_imputation", mode = "function")) {
  .imp_source_path <- file.path(.self_dir, "Imputation.R")
  if (file.exists(.imp_source_path)) {
    source(.imp_source_path, local = FALSE)
  } else {
    warning("Imputation.R not found at '", .imp_source_path,
            "'. Re-imputation will not be available.")
  }
}

# --- Benchmark methods (14 individual methods, excludes combo/softHybrid/none) ---
.IM_BENCH_METHODS <- c(
  "bpca", "knn", "mice", "missForest", "Impseq", "Impseqrob",
  "QRILC", "MLE", "MinDet", "MinProb", "PI", "min", "zero", "nbavg", "with", "limpa"
)

# =============================================================================
# SECTION 1: INTERNAL METRIC HELPERS (.im_*)
# =============================================================================

#' NRMSE — Normalized Root Mean Squared Error
#'
#' Computes NRMSE between true and imputed values at positions indicated by
#' the NA mask. Lower is better.
#'
#' @param true_mat Numeric matrix (ground truth, complete)
#' @param imp_mat  Numeric matrix (imputed from artificially-NA'd version)
#' @param na_mask  Logical matrix (TRUE = position was set to NA)
#' @return Numeric scalar (NRMSE)
#' @keywords internal
.im_nrmse <- function(true_mat, imp_mat, na_mask) {
  true_vals <- true_mat[na_mask]
  imp_vals  <- imp_mat[na_mask]
  valid <- !is.na(true_vals) & !is.na(imp_vals)
  true_vals <- true_vals[valid]
  imp_vals  <- imp_vals[valid]
  if (length(true_vals) < 2) return(NA_real_)
  v <- var(true_vals)
  if (is.na(v) || v == 0) return(NA_real_)
  sqrt(mean((imp_vals - true_vals)^2)) / sqrt(v)
}

#' Per-feature RMSE for SOR computation
#'
#' For each row (feature) with artificial NAs, computes RMSE between
#' true and imputed values. Returns a named numeric vector.
#'
#' @param true_mat Numeric matrix (ground truth)
#' @param imp_mat  Numeric matrix (imputed)
#' @param na_mask  Logical matrix
#' @return Named numeric vector of per-feature RMSE
#' @keywords internal
.im_rmse_per_feature <- function(true_mat, imp_mat, na_mask) {
  # Rows with at least one artificial NA

  row_has_na <- rowSums(na_mask) > 0
  rows_idx   <- which(row_has_na)
  if (length(rows_idx) == 0) return(numeric(0))

  rmse_vec <- vapply(rows_idx, function(i) {
    cols <- which(na_mask[i, ])
    if (length(cols) == 0) return(NA_real_)
    diffs <- imp_mat[i, cols] - true_mat[i, cols]
    valid <- !is.na(diffs)
    if (sum(valid) == 0) return(NA_real_)
    sqrt(mean(diffs[valid]^2))
  }, numeric(1))

  names(rmse_vec) <- rownames(true_mat)[rows_idx]
  rmse_vec
}

#' PSS — Procrustes Statistical Shape analysis
#'
#' PCA on both true and imputed matrices (transposed: samples as rows),
#' retain components up to 95% cumulative variance, then Procrustes SS.
#' Requires vegan (optional). Lower is better.
#'
#' @param true_mat Numeric matrix (ground truth, proteins x samples)
#' @param imp_mat  Numeric matrix (imputed)
#' @return Numeric scalar (Procrustes SS), or NA if vegan not available
#' @keywords internal
.im_pss <- function(true_mat, imp_mat) {
  if (!requireNamespace("vegan", quietly = TRUE)) return(NA_real_)
  if (ncol(true_mat) < 3 || nrow(true_mat) < 3) return(NA_real_)

  # NAguideR strategy: determine k from ground truth, apply same k to both
  pca_true <- tryCatch(
    prcomp(t(true_mat), center = TRUE, scale. = TRUE),
    error = function(e) NULL
  )
  if (is.null(pca_true)) return(NA_real_)

  cum_var <- cumsum(pca_true$sdev^2) / sum(pca_true$sdev^2)
  k <- which(cum_var >= 0.95)[1]
  if (is.na(k)) k <- length(cum_var)
  k <- max(k, 2)
  pca_true_scores <- pca_true$x[, seq_len(k), drop = FALSE]

  pca_imp <- tryCatch(
    prcomp(t(imp_mat), center = TRUE, scale. = TRUE),
    error = function(e) NULL
  )
  if (is.null(pca_imp)) return(NA_real_)

  # Use same k from ground truth (NAguideR: $x[, 1:pcazhanbi95] for both)
  k_imp <- min(k, ncol(pca_imp$x))
  if (k_imp < k) return(NA_real_)
  pca_imp_scores <- pca_imp$x[, seq_len(k), drop = FALSE]

  res <- tryCatch(
    vegan::procrustes(pca_true_scores, pca_imp_scores, symmetric = TRUE),
    error = function(e) NULL
  )

  if (is.null(res)) return(NA_real_)
  res$ss
}

#' ACC_OI — Average Correlation Coefficient of Original and Imputed
#'
#' For each feature (row) with artificial NAs, computes Pearson correlation
#' between the full true row and the full imputed row (all columns).
#' Returns the mean across features. Higher is better.
#'
#' @param true_mat Numeric matrix (ground truth)
#' @param imp_mat  Numeric matrix (imputed)
#' @param na_mask  Logical matrix
#' @return Numeric scalar (mean Pearson correlation)
#' @keywords internal
.im_acc_oi <- function(true_mat, imp_mat, na_mask) {
  row_has_na <- which(rowSums(na_mask) > 0)
  if (length(row_has_na) == 0) return(NA_real_)

  cors <- vapply(row_has_na, function(i) {
    true_row <- true_mat[i, ]
    imp_row  <- imp_mat[i, ]
    sd_true <- sd(true_row, na.rm = TRUE)
    sd_imp  <- sd(imp_row, na.rm = TRUE)
    if (is.na(sd_true) || is.na(sd_imp) || sd_true == 0 || sd_imp == 0)
      return(NA_real_)
    cor(true_row, imp_row, method = "pearson", use = "pairwise.complete.obs")
  }, numeric(1))

  mean(cors, na.rm = TRUE)
}

#' Order Method as factor (consistent with Normalization_Metrics.R pattern)
#' @param df data.frame with Method column
#' @param method_levels Character vector of method names in desired order
#' @return data.frame with Method as ordered factor
#' @keywords internal
.im_method_factor <- function(df, method_levels) {
  df$Method <- factor(df$Method, levels = method_levels)
  df
}

#' PRONE-style qualitative color palette
#' @param n Number of colors
#' @return Character vector of hex colors
#' @keywords internal
.im_prone_colors <- function(n) {
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    pal_info   <- RColorBrewer::brewer.pal.info
    qual_pals  <- pal_info[pal_info$category == "qual", ]
    col_vector <- rev(unlist(mapply(RColorBrewer::brewer.pal,
                                    qual_pals$maxcolors,
                                    rownames(qual_pals))))
    rep_len(col_vector, n)
  } else {
    grDevices::hcl(seq(15, 375, length.out = n + 1)[seq_len(n)],
                   l = 65, c = 100)
  }
}

# =============================================================================
# SECTION 2: INTRODUCTION OF ARTIFICIAL NAs
# =============================================================================

#' Introduce artificial NAs into a complete-case matrix
#'
#' Takes a matrix with no NAs (ground truth) and introduces NAs at random
#' or mimicking the NA pattern from a reference matrix.
#'
#' When `pattern = "from_data"`, replicates the NAguideR strategy
#' (Wang et al., DOI:10.1093/nar/gkz903): the proportion of affected rows
#' and the per-column NA distribution are derived entirely from `ref_mat`,
#' so `na_prop` is ignored.
#'
#' @param mat Numeric matrix (proteins x samples), no NAs (ground truth)
#' @param na_prop Numeric (0-1). Proportion of values to set as NA.
#'   Default 0.20. Only used when `pattern = "random"`;
#'   ignored when `pattern = "from_data"`.
#' @param seed Integer. Random seed. Default 42.
#' @param pattern Character: `"random"` for uniform random NAs, or
#'   `"from_data"` to replicate NAguideR's simulation strategy using
#'   the actual NA structure from `ref_mat`.
#' @param ref_mat Numeric matrix. Reference matrix with real NAs, used when
#'   `pattern = "from_data"`. Ignored otherwise.
#' @return Named list:
#'   \item{mat_with_na}{Matrix with artificial NAs introduced}
#'   \item{na_mask}{Logical matrix (TRUE = artificially set to NA)}
#'   \item{true_mat}{Original complete matrix (unchanged)}
#'   \item{summary}{List with n_total, n_na, na_rate, rows_affected,
#'     row_ratio (from_data only)}
#' @keywords internal
.im_introduce_na <- function(mat, na_prop = 0.20, seed = 42L,
                             pattern = "random", ref_mat = NULL) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  na_mask <- matrix(FALSE, nrow = nr, ncol = nc,
                    dimnames = dimnames(mat))

  if (pattern == "from_data" && !is.null(ref_mat)) {
    # --- NAguideR strategy (faithful replication) ---
    # 1. naratiox: proportion of ROWS with at least one NA in ref_mat
    incomplete   <- !complete.cases(ref_mat)
    n_incomplete <- sum(incomplete)

    if (n_incomplete == 0) {
      warning(".im_introduce_na: ref_mat has no NAs; ",
              "falling back to 'random' pattern.")
      pattern <- "random"
    } else {
      naratiox <- n_incomplete / nrow(ref_mat)

      # 2. nacolratio: per-column NA rate among incomplete rows
      ref_incomplete <- ref_mat[incomplete, , drop = FALSE]
      nacolratio <- colSums(is.na(ref_incomplete)) / n_incomplete

      # 3. Select which rows in ground-truth will receive NAs
      nanum <- round(naratiox * nr)
      nanum <- max(1L, min(nanum, nr - 1L))  # at least 1, keep at least 1 clean
      set.seed(seed)
      samplenaindex <- sample.int(nr, nanum)

      # 4. Per column: sample a subset of those rows and set to NA
      for (j in seq_len(nc)) {
        n_na_col <- round(nacolratio[j] * nanum)
        if (n_na_col > 0 && n_na_col <= length(samplenaindex)) {
          set.seed(seed + j)  # per-column reproducibility (NAguideR: set.seed(i))
          rows_j <- sample(samplenaindex, n_na_col)
          na_mask[rows_j, j] <- TRUE
        }
      }

      # 5. Fallback: rows selected but received no NAs across all columns
      #    (can happen if nacolratio is very low). NAguideR assigns NAs
      #    row-wise in that case using naratiox * ncol as per-row count.
      rows_no_na <- samplenaindex[rowSums(na_mask[samplenaindex, , drop = FALSE]) == 0]
      if (length(rows_no_na) > 0) {
        eachrownaratio <- max(1L, round(naratiox * nc))
        for (idx in seq_along(rows_no_na)) {
          set.seed(seed + nc + idx)
          cols_k <- sample.int(nc, min(eachrownaratio, nc))
          na_mask[rows_no_na[idx], cols_k] <- TRUE
        }
      }
    }
  }

  if (pattern == "random") {
    na_prop <- na_prop %||% 0.20
    set.seed(seed)
    n_total <- nr * nc
    n_na    <- round(n_total * na_prop)
    idx     <- sample.int(n_total, n_na)
    na_mask[idx] <- TRUE
  }

  mat_with_na <- mat
  mat_with_na[na_mask] <- NA_real_

  summary_list <- list(
    n_total = nr * nc,
    n_na    = sum(na_mask),
    na_rate = mean(na_mask)
  )
  if (pattern == "from_data" && exists("naratiox", inherits = FALSE)) {
    summary_list$rows_affected <- nanum
    summary_list$row_ratio     <- naratiox
  }

  list(
    mat_with_na = mat_with_na,
    na_mask     = na_mask,
    true_mat    = mat,
    summary     = summary_list
  )
}

# =============================================================================
# SECTION 3: RE-IMPUTATION ENGINE
# =============================================================================

#' Re-impute a matrix with multiple methods
#'
#' Iterates over a vector of imputation method names, calling
#' `.dispatch_imputation()` from Imputation.R for each. Also supports
#' combo and softHybrid methods via `combo_methods` parameter.
#' Methods that fail are omitted with a warning.
#'
#' @param mat_with_na Numeric matrix with artificial NAs
#' @param methods Character vector of individual method names
#' @param combo_methods Named list of combo/softHybrid configurations. Each
#'   element is a list with: `mode` ("combo" or "softHybrid", default "combo"),
#'   `mar_method`, `mnar_method`, and optionally other parameters. The list
#'   name is used as the method label.
#' @param condition Factor or character vector of conditions (required for
#'   combo methods, ignored for individual methods).
#' @param method_args Named list of per-method argument lists
#' @param with_value Constant for method "with"
#' @param verbose Logical. Print progress. Default TRUE.
#' @return Named list of imputed matrices (one per successful method)
#' @keywords internal
.im_reimpute <- function(mat_with_na, methods = character(0),
                         combo_methods = list(),
                         condition = NULL,
                         method_args = list(),
                         with_value = NA_real_, verbose = TRUE) {
  if (!exists(".dispatch_imputation", mode = "function")) {
    stop(".dispatch_imputation() not available. ",
         "Ensure Imputation.R is sourced before using this function.")
  }

  results <- list()

  # --- Individual methods ---
  for (m in methods) {
    if (verbose) message("  Re-imputing with method: ", m, " ...")
    results[[m]] <- tryCatch(
      .dispatch_imputation(mat_with_na, m, method_args, with_value),
      error = function(e) {
        warning("Method '", m, "' failed: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
  }

  # --- Combo / softHybrid methods ---
  for (label in names(combo_methods)) {
    cfg  <- combo_methods[[label]]
    mode <- cfg$mode %||% "combo"
    mar  <- cfg$mar_method %||% "Impseqrob"
    mnar <- cfg$mnar_method %||% "min"

    if (verbose) message("  Re-imputing with ", mode, " (", mar, "+", mnar, ") as '", label, "' ...")

    results[[label]] <- tryCatch({
      if (mode == "combo") {
        if (is.null(condition))
          stop("condition is required for combo methods.")
        if (!exists(".impute_combo", mode = "function"))
          stop(".impute_combo() not available. Ensure Imputation.R is sourced.")
        res <- .impute_combo(
          x           = mat_with_na,
          condition   = condition,
          mar_method  = mar,
          mnar_method = mnar,
          method_args = method_args,
          with_value  = with_value
        )
        res$x_imputed
      } else if (mode == "softHybrid") {
        if (!exists(".impute_softHybrid", mode = "function"))
          stop(".impute_softHybrid() not available. Ensure Imputation.R is sourced.")
        sh_args <- cfg[setdiff(names(cfg), c("mode", "mar_method", "mnar_method"))]
        res <- do.call(.impute_softHybrid, c(
          list(x = mat_with_na, mar_method = mar, mnar_method = mnar,
               method_args = method_args, with_value = with_value),
          sh_args
        ))
        res$x_imputed
      } else {
        stop("Unknown mode '", mode, "'. Use 'combo' or 'softHybrid'.")
      }
    }, error = function(e) {
      warning("Combo method '", label, "' failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
  }

  # Remove failed methods
  results[!vapply(results, is.null, logical(1))]
}

# =============================================================================
# SECTION 4: IMPORT FUNCTION
# =============================================================================

#' Import imputed matrices into a SummarizedExperiment
#'
#' Reads all TSV files matching a pattern from a directory, one assay per file.
#' Each TSV must have column 1 as protein identifier and remaining columns as
#' sample names matching the `sample_col` in the metadata file.
#'
#' @param tsv_dir Character. Directory containing the TSV files.
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
#' se_imp <- import_imp_matrices(
#'   tsv_dir       = "./results",
#'   metadata_path = "./data/metadata.tsv"
#' )
#' SummarizedExperiment::assayNames(se_imp)
#' }
#' @export
import_imp_matrices <- function(tsv_dir,
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

  # --- Helper: read TSV ---
  .read_tsv_file <- function(path) {
    if (requireNamespace("readr", quietly = TRUE)) {
      as.data.frame(readr::read_tsv(path, show_col_types = FALSE),
                    stringsAsFactors = FALSE)
    } else {
      read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    }
  }

  # --- Read all files to get proteins and sample columns ---
  all_dfs <- vector("list", length(files))
  for (i in seq_along(files)) {
    all_dfs[[i]] <- .read_tsv_file(files[i])
  }

  protein_col <- colnames(all_dfs[[1]])[1]
  sample_cols <- setdiff(colnames(all_dfs[[1]]), protein_col)

  # --- Validate samples against metadata ---
  missing_meta <- setdiff(sample_cols, meta[[sample_col]])
  if (length(missing_meta) > 0)
    stop("Samples in TSV not found in metadata: ",
         paste(missing_meta, collapse = ", "))

  meta_aligned <- meta[sample_cols, , drop = FALSE]
  rownames(meta_aligned) <- sample_cols

  # --- Find common proteins across all files ---
  protein_sets <- lapply(all_dfs, function(df) df[[protein_col]])
  proteins <- Reduce(intersect, protein_sets)

  if (length(proteins) == 0)
    stop("No common proteins found across all TSV files.")

  n_orig <- vapply(protein_sets, length, integer(1))
  if (any(n_orig != length(proteins))) {
    message("import_imp_matrices: files have different row counts (",
            paste(unique(n_orig), collapse = ", "),
            "). Using intersection: ", length(proteins), " common proteins.")
  }

  # --- Build assay list (aligned to common proteins) ---
  assay_list <- vector("list", length(files))
  names(assay_list) <- method_names

  for (i in seq_along(files)) {
    df  <- all_dfs[[i]]
    rownames(df) <- df[[protein_col]]
    df_sub <- df[proteins, sample_cols, drop = FALSE]
    mat <- as.matrix(df_sub)
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

  message("import_imp_matrices: loaded ", length(files), " assay(s) — ",
          paste(method_names, collapse = ", "))
  message("  Proteins: ", nrow(se), " | Samples: ", ncol(se))
  se
}

# =============================================================================
# SECTION 4b: PREPARE SE FROM PREPROCESSING + BEST NORMALIZATION
# =============================================================================

#' Build a SummarizedExperiment for imputation benchmarking
#'
#' Convenience function that creates a SE from a `spectronaut_data` object and
#' applies the winning normalization method (from normalization benchmarking).
#' The result is ready for `imputation_metrics()`.
#'
#' Internally calls `nm_prepare_se()` to build the baseline SE with assay
#' `"log2"`, then applies the chosen normalization via
#' `.nm_dispatch_normalization()`.
#'
#' @param preprocessing `spectronaut_data` list from `preprocess_spectronaut()`.
#' @param norm_method Character scalar. Normalization method name (e.g.
#'   `"cycloess"`). If provided, used directly. Default `NULL`.
#' @param pc1_rank data.frame from `nm_rank_pc1()` or
#'   `normalization_metrics()$pc1_rank`. The first row's `Method` column is
#'   used as winner. Ignored if `norm_method` is provided. Default `NULL`.
#' @param min_reps Minimum replicates with non-NA values per group for protein
#'   filtering. If NULL, auto-computed as half the smallest group. Default `NULL`.
#' @param min_groups Minimum groups meeting `min_reps` (default: 1).
#' @param covariate_df Optional covariate data.frame for paired designs
#'   (must contain a `Column` column). Default `NULL`.
#' @param norm_method_args Named list of per-method arguments for normalization.
#'   Default `list()`.
#' @param include_log2 Logical. Include the baseline `"log2"` assay in the
#'   returned SE. Default `TRUE`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @return SummarizedExperiment with assays `"log2"` (optional) + winner method.
#'
#' @examples
#' \dontrun{
#' source("R/Imputation_Metrics.R")
#'
#' # ---- Option A: known method ----
#' se_imp <- im_prepare_se(preprocessing, norm_method = "cycloess")
#'
#' # ---- Option B: auto-pick from pc1_rank ----
#' # (after running normalization_metrics() in Normalization_Metrics.R)
#' se_imp <- im_prepare_se(preprocessing, pc1_rank = nm_res$pc1_rank)
#'
#' SummarizedExperiment::assayNames(se_imp)  # "log2", "cycloess"
#' res <- imputation_metrics(se_imp, assay_name = "cycloess")
#' }
#' @export
im_prepare_se <- function(preprocessing,
                          norm_method      = NULL,
                          pc1_rank         = NULL,
                          min_reps         = NULL,
                          min_groups       = 1,
                          covariate_df     = NULL,
                          norm_method_args = list(),
                          include_log2     = TRUE,
                          verbose          = TRUE) {

  # --- Validate preprocessing ---
  if (!inherits(preprocessing, "spectronaut_data"))
    stop("'preprocessing' must be a spectronaut_data object ",
         "(output of preprocess_spectronaut()).")

  # --- Resolve winner method ---
  if (is.null(norm_method) && is.null(pc1_rank))
    stop("At least one of 'norm_method' or 'pc1_rank' must be provided.")

  if (is.null(norm_method)) {
    if (!is.data.frame(pc1_rank) || !"Method" %in% colnames(pc1_rank) ||
        nrow(pc1_rank) == 0)
      stop("'pc1_rank' must be a data.frame with a non-empty 'Method' column ",
           "(output of nm_rank_pc1() or normalization_metrics()$pc1_rank).")
    norm_method <- as.character(pc1_rank$Method[1])
    if (verbose)
      message("im_prepare_se: winner from pc1_rank -> '", norm_method, "'")
  }

  # --- Lazy-source Normalization_Metrics.R ---
  if (!exists("nm_prepare_se", mode = "function") ||
      !exists(".nm_dispatch_normalization", mode = "function")) {
    nm_source_path <- file.path(.self_dir, "Normalization_Metrics.R")
    if (file.exists(nm_source_path)) {
      source(nm_source_path, local = FALSE)
    } else {
      stop("Normalization_Metrics.R not found at '", nm_source_path,
           "'. Required for im_prepare_se().")
    }
  }

  # --- Build baseline SE with assay "log2" ---
  if (verbose) message("im_prepare_se: building baseline SE ...")
  se <- nm_prepare_se(preprocessing,
                      min_reps     = min_reps,
                      min_groups   = min_groups,
                      covariate_df = covariate_df,
                      verbose      = verbose)

  # --- Apply winner normalization ---
  if (verbose) message("im_prepare_se: applying normalization '", norm_method, "' ...")
  x_log2 <- SummarizedExperiment::assay(se, "log2")
  x_norm <- .nm_dispatch_normalization(x_log2, norm_method, norm_method_args)

  # --- Build final SE ---
  if (include_log2) {
    assay_list <- list(x_log2, x_norm)
    names(assay_list) <- c("log2", norm_method)
  } else {
    assay_list <- list(x_norm)
    names(assay_list) <- norm_method
  }

  se_out <- SummarizedExperiment::SummarizedExperiment(
    assays   = assay_list,
    colData  = SummarizedExperiment::colData(se),
    rowData  = SummarizedExperiment::rowData(se),
    metadata = S4Vectors::metadata(se)
  )

  if (verbose) {
    message("im_prepare_se: done. Assays: ",
            paste(SummarizedExperiment::assayNames(se_out), collapse = ", "))
    message("  Suggestion: imputation_metrics(se, assay_name = \"",
            norm_method, "\")")
  }

  se_out
}

# =============================================================================
# SECTION 5: COMPUTE METRICS
# =============================================================================

#' Compute imputation quality metrics (NRMSE, SOR, PSS, ACC_OI)
#'
#' Extracts complete-case rows from the assay, introduces artificial NAs,
#' re-imputes with each method, and computes four quality metrics plus
#' a combined ranking.
#'
#' @param se SummarizedExperiment (from `import_imp_matrices()` or pipeline)
#' @param assay_name Character. Name of the assay to use as starting point
#'   (typically the normalized assay, pre-imputation, with real NAs).
#' @param methods Character vector of individual imputation methods to benchmark.
#'   Default: `.IM_BENCH_METHODS` (14 methods). Use `character(0)` or `NULL`
#'   to skip individual methods.
#' @param combo_methods Named list of combo/softHybrid configurations. Each
#'   element is a list with: `mode` ("combo" or "softHybrid", default "combo"),
#'   `mar_method`, `mnar_method`. The list name is used as the method label
#'   in results. Default: empty list (no combo methods).
#' @param condition_col Character. Column in colData with condition labels,
#'   required for combo methods. Default `"Condition"`.
#' @param na_prop Numeric (0-1). Proportion of artificial NAs. Default 0.20.
#'   Only used when `pattern = "random"`. Ignored when `pattern = "from_data"`
#'   (NA proportion is derived from the data, replicating NAguideR).
#' @param seed Integer. Random seed. Default 42.
#' @param pattern Character. NA introduction pattern: `"random"` for uniform
#'   random NAs, or `"from_data"` to replicate NAguideR's simulation strategy
#'   (row ratio + per-column distribution from the actual data).
#' @param method_args Named list of per-method argument lists.
#' @param with_value Constant value for method `"with"`.
#' @param verbose Logical. Print progress. Default TRUE.
#' @return data.frame with columns: Method, NRMSE, SOR, PSS, ACC_OI,
#'   NRMSE_Rank, SOR_Rank, PSS_Rank, ACC_OI_Rank, Rank_Mean.
#'   Ordered by Rank_Mean (best first).
#'
#' @examples
#' \dontrun{
#' # Individual methods only
#' metrics <- im_compute_metrics(se, assay_name = "cycloess",
#'                                methods = c("knn", "min", "zero"))
#'
#' # Combo + individual
#' metrics <- im_compute_metrics(se, assay_name = "cycloess",
#'   methods = c("knn", "min"),
#'   combo_methods = list(
#'     "Impseq+min" = list(mar_method = "Impseq", mnar_method = "min")
#'   ))
#' }
#' @export
im_compute_metrics <- function(se,
                               assay_name    = NULL,
                               methods       = .IM_BENCH_METHODS,
                               combo_methods = list(),
                               condition_col = "Condition",
                               na_prop       = 0.20,
                               seed          = 42L,
                               pattern       = "random",
                               method_args   = list(),
                               with_value    = NA_real_,
                               verbose       = TRUE) {
  # --- Validate ---
  stopifnot(inherits(se, "SummarizedExperiment"))
  all_assays <- SummarizedExperiment::assayNames(se)

  if (is.null(assay_name)) {
    assay_name <- all_assays[1]
    if (verbose) message("im_compute_metrics: using first assay '", assay_name, "'")
  }
  if (!assay_name %in% all_assays)
    stop("Assay '", assay_name, "' not found in SE. Available: ",
         paste(all_assays, collapse = ", "))

  # Inform about optional packages
  if (!requireNamespace("vegan", quietly = TRUE))
    message("im_compute_metrics: 'vegan' not installed — PSS will be NA.")

  # --- Extract matrix and complete cases ---
  mat_full <- SummarizedExperiment::assay(se, assay_name)
  complete_rows <- complete.cases(mat_full)
  mat_complete  <- mat_full[complete_rows, , drop = FALSE]

  if (nrow(mat_complete) < 10) {
    stop("Only ", nrow(mat_complete), " complete-case rows available. ",
         "Need at least 10 for reliable metric computation.")
  }

  if (verbose) {
    message("im_compute_metrics: ", nrow(mat_complete), " complete-case rows out of ",
            nrow(mat_full), " total.")
  }

  # --- Introduce artificial NAs ---
  ref_mat <- if (pattern == "from_data") mat_full else NULL
  na_result <- .im_introduce_na(mat_complete, na_prop = na_prop, seed = seed,
                                pattern = pattern, ref_mat = ref_mat)

  if (verbose) {
    message("  Artificial NAs introduced: ", na_result$summary$n_na,
            " (", round(na_result$summary$na_rate * 100, 1), "%)")
    if (pattern == "from_data" && !is.null(na_result$summary$row_ratio)) {
      message("  NAguideR strategy: ", na_result$summary$rows_affected,
              " rows affected (",
              round(na_result$summary$row_ratio * 100, 1),
              "% row ratio from data)")
    }
  }

  # --- Resolve condition vector (needed for combo methods) ---
  condition <- NULL
  if (length(combo_methods) > 0) {
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    if (!condition_col %in% colnames(cd))
      stop("Column '", condition_col, "' not found in colData. ",
           "Required for combo methods.")
    condition <- as.factor(cd[[condition_col]])
  }

  # --- Exclude "none" (it doesn't impute, metrics are meaningless) ---
  methods <- methods %||% character(0)
  if ("none" %in% methods) {
    warning("Method 'none' excluded from imputation metrics ",
            "(it does not impute; metrics are not applicable).",
            call. = FALSE)
    methods <- setdiff(methods, "none")
  }
  imp_results <- .im_reimpute(na_result$mat_with_na,
                              methods       = methods,
                              combo_methods = combo_methods,
                              condition     = condition,
                              method_args   = method_args,
                              with_value    = with_value,
                              verbose       = verbose)

  if (length(imp_results) == 0)
    stop("All imputation methods failed. Cannot compute metrics.")

  successful_methods <- names(imp_results)

  # --- Compute per-method metrics ---
  true_mat <- na_result$true_mat
  na_mask  <- na_result$na_mask

  # Step 1: Per-feature RMSE for all methods (for SOR)
  rmse_list <- lapply(imp_results, function(imp_mat) {
    .im_rmse_per_feature(true_mat, imp_mat, na_mask)
  })

  # Step 2: Build feature x method RMSE matrix for SOR
  all_features <- unique(unlist(lapply(rmse_list, names)))
  rmse_matrix  <- matrix(NA_real_, nrow = length(all_features),
                         ncol = length(successful_methods),
                         dimnames = list(all_features, successful_methods))
  for (m in successful_methods) {
    rmse_vec <- rmse_list[[m]]
    common   <- intersect(names(rmse_vec), all_features)
    rmse_matrix[common, m] <- rmse_vec[common]
  }

  # Step 3: Rank per feature, then sum ranks → SOR
  if (length(successful_methods) == 1L) {
    # Single method: all ranks = 1, SOR = number of features
    sor_vec <- setNames(nrow(rmse_matrix), successful_methods)
  } else {
    rank_per_feature <- t(apply(rmse_matrix, 1, function(row) {
      rank(row, na.last = "keep", ties.method = "average")
    }))
    sor_vec <- colSums(rank_per_feature, na.rm = TRUE)
  }

  # Step 4: All four metrics
  metrics_rows <- vector("list", length(successful_methods))
  for (i in seq_along(successful_methods)) {
    m <- successful_methods[i]
    imp_mat <- imp_results[[m]]
    metrics_rows[[i]] <- data.frame(
      Method = m,
      NRMSE  = .im_nrmse(true_mat, imp_mat, na_mask),
      SOR    = unname(sor_vec[m]),
      PSS    = .im_pss(true_mat, imp_mat),
      ACC_OI = .im_acc_oi(true_mat, imp_mat, na_mask),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  metrics_df <- do.call(rbind, metrics_rows)

  # --- Compute ranks ---
  n <- nrow(metrics_df)

  # Lower is better: NRMSE, SOR, PSS
  metrics_df$NRMSE_Rank <- rank(metrics_df$NRMSE, na.last = "keep",
                                ties.method = "average")
  metrics_df$SOR_Rank   <- rank(metrics_df$SOR, na.last = "keep",
                                ties.method = "average")
  metrics_df$PSS_Rank   <- rank(metrics_df$PSS, na.last = "keep",
                                ties.method = "average")

  # Higher is better: ACC_OI → rank descending
  if (all(is.na(metrics_df$ACC_OI))) {
    metrics_df$ACC_OI_Rank <- rep(NA_real_, n)
  } else {
    metrics_df$ACC_OI_Rank <- rank(-metrics_df$ACC_OI, na.last = "keep",
                                   ties.method = "average")
  }

  # --- Rank_Mean ---
  rank_cols <- c("NRMSE_Rank", "SOR_Rank", "PSS_Rank", "ACC_OI_Rank")
  metrics_df$Rank_Mean <- rowMeans(metrics_df[, rank_cols], na.rm = TRUE)

  # --- Order by Rank_Mean ---
  metrics_df <- metrics_df[order(metrics_df$Rank_Mean), ]
  rownames(metrics_df) <- NULL

  metrics_df
}

# =============================================================================
# SECTION 6: PLOT FUNCTIONS (6 plots, ggplot2)
# =============================================================================

# --------------------------------------------------------------------------
# 1. NRMSE bar plot
# --------------------------------------------------------------------------

#' NRMSE bar plot per imputation method
#'
#' Bar chart of NRMSE values, ordered ascending (lower is better).
#'
#' @param metrics_df data.frame from `im_compute_metrics()`.
#' @param ... Additional arguments (unused).
#' @return ggplot object.
#' @export
im_plot_nrmse <- function(metrics_df, ...) {
  df <- metrics_df[order(metrics_df$NRMSE), ]
  df$Method <- factor(df$Method, levels = df$Method)
  col_vector <- .im_prone_colors(nrow(df))

  ggplot2::ggplot(df, ggplot2::aes(x = Method, y = NRMSE, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_label(
      ggplot2::aes(label = sprintf("%.3f", NRMSE)),
      size = 2.8, fill = "white", linewidth = 0.2, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = col_vector) +
    ggplot2::labs(
      title = "NRMSE per Imputation Method (lower is better)",
      x = NULL, y = "NRMSE"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      plot.title  = ggplot2::element_text(face = "bold", size = 12)
    )
}

# --------------------------------------------------------------------------
# 2. SOR bar plot
# --------------------------------------------------------------------------

#' SOR bar plot per imputation method
#'
#' Bar chart of Sum of Ranks (SOR), ordered ascending (lower is better).
#'
#' @inheritParams im_plot_nrmse
#' @return ggplot object.
#' @export
im_plot_sor <- function(metrics_df, ...) {
  df <- metrics_df[order(metrics_df$SOR), ]
  df$Method <- factor(df$Method, levels = df$Method)
  col_vector <- .im_prone_colors(nrow(df))

  ggplot2::ggplot(df, ggplot2::aes(x = Method, y = SOR, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_label(
      ggplot2::aes(label = sprintf("%.0f", SOR)),
      size = 2.8, fill = "white", linewidth = 0.2, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = col_vector) +
    ggplot2::labs(
      title = "SOR per Imputation Method (lower is better)",
      x = NULL, y = "Sum of Ranks (SOR)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      plot.title  = ggplot2::element_text(face = "bold", size = 12)
    )
}

# --------------------------------------------------------------------------
# 3. PSS bar plot
# --------------------------------------------------------------------------

#' PSS bar plot per imputation method
#'
#' Bar chart of Procrustes Statistical Shape (PSS), ordered ascending
#' (lower is better). Requires vegan; NA values shown as 0-height bars.
#'
#' @inheritParams im_plot_nrmse
#' @return ggplot object.
#' @export
im_plot_pss <- function(metrics_df, ...) {
  df <- metrics_df[order(metrics_df$PSS), ]
  df$Method <- factor(df$Method, levels = df$Method)
  col_vector <- .im_prone_colors(nrow(df))

  ggplot2::ggplot(df, ggplot2::aes(x = Method, y = PSS, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_label(
      ggplot2::aes(label = ifelse(is.na(PSS), "NA", sprintf("%.4f", PSS))),
      size = 2.8, fill = "white", linewidth = 0.2, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = col_vector) +
    ggplot2::labs(
      title = "PSS per Imputation Method (lower is better)",
      subtitle = "Procrustes Statistical Shape",
      x = NULL, y = "PSS"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      plot.title  = ggplot2::element_text(face = "bold", size = 12)
    )
}

# --------------------------------------------------------------------------
# 4. ACC_OI bar plot
# --------------------------------------------------------------------------

#' ACC_OI bar plot per imputation method
#'
#' Bar chart of Average Correlation Coefficient (ACC_OI), ordered descending
#' (higher is better).
#'
#' @inheritParams im_plot_nrmse
#' @return ggplot object.
#' @export
im_plot_acc_oi <- function(metrics_df, ...) {
  df <- metrics_df[order(-metrics_df$ACC_OI), ]
  df$Method <- factor(df$Method, levels = df$Method)
  col_vector <- .im_prone_colors(nrow(df))

  ggplot2::ggplot(df, ggplot2::aes(x = Method, y = ACC_OI, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_label(
      ggplot2::aes(label = sprintf("%.3f", ACC_OI)),
      size = 2.8, fill = "white", linewidth = 0.2, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = col_vector) +
    ggplot2::labs(
      title = "ACC_OI per Imputation Method (higher is better)",
      subtitle = "Average Correlation Coefficient (Original vs Imputed)",
      x = NULL, y = "ACC_OI"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      plot.title  = ggplot2::element_text(face = "bold", size = 12)
    )
}

# --------------------------------------------------------------------------
# 5. Ranking heatmap
# --------------------------------------------------------------------------

#' Ranking heatmap of imputation methods
#'
#' Heatmap (geom_tile + geom_text) showing per-metric ranks and the
#' combined Rank_Mean. Rows = methods (ordered by Rank_Mean), columns = metrics.
#'
#' @inheritParams im_plot_nrmse
#' @return ggplot object.
#' @export
im_plot_ranking <- function(metrics_df, ...) {
  # Order by Rank_Mean
  df <- metrics_df[order(metrics_df$Rank_Mean), ]
  method_order <- df$Method

  # Pivot rank columns to long format
  rank_cols <- c("NRMSE_Rank", "SOR_Rank", "PSS_Rank", "ACC_OI_Rank", "Rank_Mean")
  rank_df <- df[, c("Method", rank_cols)]

  long_df <- tidyr::pivot_longer(
    rank_df,
    cols      = tidyr::all_of(rank_cols),
    names_to  = "Metric",
    values_to = "Rank"
  )

  long_df$Method <- factor(long_df$Method, levels = rev(method_order))
  long_df$Metric <- factor(long_df$Metric, levels = rank_cols)

  ggplot2::ggplot(long_df,
    ggplot2::aes(x = Metric, y = Method, fill = Rank)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(Rank), "NA", sprintf("%.1f", Rank))),
      size = 3.5, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = median(long_df$Rank, na.rm = TRUE),
      na.value = "grey80",
      name = "Rank") +
    ggplot2::labs(
      title = "Imputation Method Ranking",
      subtitle = "Lower rank = better performance",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y  = ggplot2::element_text(size = 9),
      plot.title   = ggplot2::element_text(face = "bold", size = 12),
      panel.grid   = ggplot2::element_blank()
    )
}

# --------------------------------------------------------------------------
# 6. Faceted metrics bar chart
# --------------------------------------------------------------------------

#' Faceted bar chart of all 4 imputation quality metrics
#'
#' Similar to `nm_plot_metrics()`: one facet per metric with free y-scales,
#' labeled values on bars.
#'
#' @inheritParams im_plot_nrmse
#' @return ggplot object.
#' @export
im_plot_metrics <- function(metrics_df, ...) {
  value_cols <- c("NRMSE", "SOR", "PSS", "ACC_OI")
  present_cols <- intersect(value_cols, colnames(metrics_df))

  long_df <- tidyr::pivot_longer(
    metrics_df[, c("Method", present_cols)],
    cols      = tidyr::all_of(present_cols),
    names_to  = "Metric",
    values_to = "Value"
  )

  # Order methods by Rank_Mean
  method_order <- metrics_df$Method[order(metrics_df$Rank_Mean)]
  long_df$Method <- factor(long_df$Method, levels = method_order)
  long_df$Metric <- factor(long_df$Metric, levels = value_cols)
  col_vector <- .im_prone_colors(length(method_order))

  ggplot2::ggplot(long_df,
    ggplot2::aes(x = Method, y = Value, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_label(
      ggplot2::aes(label = ifelse(is.na(Value), "NA",
                                  sprintf("%.3f", Value))),
      size = 2.2, fill = "white", linewidth = 0.2, na.rm = TRUE) +
    ggplot2::facet_wrap(~ Metric, ncol = 2, scales = "free_y") +
    ggplot2::scale_fill_manual(values = col_vector) +
    ggplot2::labs(
      title = "Imputation Quality Metrics per Method",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text  = ggplot2::element_text(face = "bold")
    )
}

# =============================================================================
# SECTION 7: MAIN ORCHESTRATOR
# =============================================================================

#' Generate imputation quality metric plots and ranking
#'
#' Computes NRMSE, SOR, PSS, ACC_OI via ground-truth simulation, then
#' generates up to 6 diagnostic plots. Each plot is wrapped in `tryCatch()`
#' so that a failure in one does not abort all.
#'
#' @param se SummarizedExperiment (from `import_imp_matrices()` or pipeline).
#' @param assay_name Character. Assay name to use as starting point.
#'   NULL (default) = first assay.
#' @param methods Character vector of individual methods to benchmark.
#'   Default: `.IM_BENCH_METHODS` (14 methods). Use `NULL` or `character(0)`
#'   to skip individual methods when only using combo_methods.
#' @param combo_methods Named list of combo/softHybrid configurations.
#'   Each element: `list(mar_method, mnar_method, mode)`.
#'   See `im_compute_metrics()` for details.
#' @param condition_col Character. Column in colData for combo methods.
#'   Default `"Condition"`.
#' @param na_prop Numeric (0-1). Proportion of artificial NAs. Default 0.20.
#'   Only used when `pattern = "random"`. Ignored when `pattern = "from_data"`
#'   (NA proportion is derived from the data, replicating NAguideR).
#' @param seed Integer. Random seed. Default 42.
#' @param pattern Character. `"random"` or `"from_data"`. Default `"random"`.
#' @param method_args Named list of per-method argument lists.
#' @param with_value Constant for method `"with"`.
#' @param plots Character vector of plot names or `"all"` (default).
#'   Valid: `"nrmse"`, `"sor"`, `"pss"`, `"acc_oi"`, `"ranking"`, `"metrics"`.
#' @param verbose Logical. Print progress. Default TRUE.
#' @return Named list of ggplot objects (or NULL for failed plots), plus
#'   `metrics_table`: a data.frame from `im_compute_metrics()`.
#'
#' @examples
#' \dontrun{
#' # ---- 1. From pipeline SE ----
#' res <- imputation_metrics(se, assay_name = "cycloess",
#'   methods = c("knn", "min"),
#'   combo_methods = list(
#'     "Impseq+min" = list(mar_method = "Impseq", mnar_method = "min")
#'   ))
#' res$metrics_table
#' res$ranking
#'
#' # ---- 1b. From preprocessing + best normalization (via im_prepare_se) ----
#'
#' # Option A: known method
#' se_imp <- im_prepare_se(preprocessing, norm_method = "cycloess")
#'
#' # Option B: auto-pick from pc1_rank
#' # (after running normalization_metrics() in Normalization_Metrics.R)
#' se_imp <- im_prepare_se(preprocessing, pc1_rank = nm_res$pc1_rank)
#'
#' SummarizedExperiment::assayNames(se_imp)  # "log2", "cycloess"
#' res <- imputation_metrics(se_imp, assay_name = "cycloess")
#' }
#' @export
imputation_metrics <- function(se,
                               assay_name    = NULL,
                               methods       = .IM_BENCH_METHODS,
                               combo_methods = list(),
                               condition_col = "Condition",
                               na_prop       = 0.20,
                               seed          = 42L,
                               pattern       = "random",
                               method_args   = list(),
                               with_value    = NA_real_,
                               plots         = "all",
                               verbose       = TRUE) {
  # --- Required packages check ---
  for (pkg in c("ggplot2", "dplyr", "tidyr", "SummarizedExperiment", "S4Vectors")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for imputation_metrics().")
  }

  # --- Compute metrics ---
  if (verbose) message("=== COMPUTING IMPUTATION METRICS ===")
  metrics_df <- im_compute_metrics(
    se            = se,
    assay_name    = assay_name,
    methods       = methods,
    combo_methods = combo_methods,
    condition_col = condition_col,
    na_prop       = na_prop,
    seed          = seed,
    pattern       = pattern,
    method_args   = method_args,
    with_value    = with_value,
    verbose       = verbose
  )

  # --- Plot registry ---
  all_plot_names <- c("nrmse", "sor", "pss", "acc_oi", "ranking", "metrics")

  plot_fns <- list(
    nrmse   = function() im_plot_nrmse(metrics_df),
    sor     = function() im_plot_sor(metrics_df),
    pss     = function() im_plot_pss(metrics_df),
    acc_oi  = function() im_plot_acc_oi(metrics_df),
    ranking = function() im_plot_ranking(metrics_df),
    metrics = function() im_plot_metrics(metrics_df)
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
        warning("im_plot_", nm, "() failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  # --- Always include metrics_table ---
  result[["metrics_table"]] <- metrics_df

  n_ok   <- sum(!sapply(result[setdiff(names(result), "metrics_table")], is.null))
  n_fail <- length(selected) - n_ok
  if (verbose) {
    message("imputation_metrics: ", n_ok, " plot(s) generated",
            if (n_fail > 0) paste0(", ", n_fail, " failed") else ".")
  }

  result
}

# =============================================================================
# SECTION 8: EXAMPLE WORKFLOW
# =============================================================================
#
# Assumes:
#   - ./results/ contains normalized matrix TSV files (matrix_log2_*.tsv)
#   - ./data/metadata.tsv has at least: Column, Condition
#   - Imputation.R is accessible in ./R/
#
# Run with:
#   source("R/Imputation_Metrics.R")
# -----------------------------------------------------------------------------

if (FALSE) {

  # ---- 1. Load normalized matrices into a SummarizedExperiment ---------------

  se_nm <- import_norm_matrices(
    tsv_dir       = "./results",
    metadata_path = "./data/metadata.tsv",
    pattern       = "matrix_log2_.*\\.tsv$"
  )

  SummarizedExperiment::assayNames(se_nm)
  dim(se_nm)


  # ---- 2. Generate all 6 plots + metrics table (full benchmark) --------------

  res <- imputation_metrics(se_nm, assay_name = "cycloess")

  names(res)  # nrmse sor pss acc_oi ranking metrics metrics_table


  # ---- 3. Inspect individual plots -------------------------------------------

  res$nrmse           # NRMSE bar chart (lower = better)
  res$sor             # SOR bar chart (lower = better)
  res$pss             # PSS bar chart (lower = better, requires vegan)
  res$acc_oi          # ACC_OI bar chart (higher = better)
  res$ranking         # Heatmap of ranks per metric
  res$metrics         # Faceted bar chart of all 4 metrics

  # Metrics table (data.frame, always present)
  res$metrics_table


  # ---- 4. Subset of methods (faster) -----------------------------------------

  res_fast <- imputation_metrics(
    se_nm,
    assay_name = "cycloess",
    methods    = c("knn", "Impseqrob", "QRILC", "min", "zero")
  )
  res_fast$metrics_table
  res_fast$ranking


  # ---- 5. Standalone metrics computation (no plots) --------------------------

  metrics_df <- im_compute_metrics(
    se_nm,
    assay_name = "cycloess",
    methods    = c("knn", "min", "MinDet", "zero"),
    verbose    = TRUE
  )
  print(metrics_df)


  # ---- 6. Individual plot functions ------------------------------------------

  im_plot_nrmse(metrics_df)
  im_plot_ranking(metrics_df)


  # ---- 7. Combo and softHybrid methods ---------------------------------------

  # Compare individual methods alongside combo (MAR+MNAR) strategies
  res_combo <- imputation_metrics(
    se_nm,
    assay_name = "cycloess",
    methods    = c("knn", "min", "MinDet"),
    combo_methods = list(
      "Impseq+min"      = list(mar_method = "Impseq", mnar_method = "min"),
      "Impseqrob+min"   = list(mar_method = "Impseqrob", mnar_method = "min"),
      "bpca+MinProb"    = list(mar_method = "bpca", mnar_method = "MinProb"),
      "sH_knn+QRILC"    = list(mode = "softHybrid",
                                mar_method = "knn", mnar_method = "QRILC")
    )
  )
  res_combo$metrics_table
  res_combo$ranking

  # Only combo methods (no individual)
  res_combo_only <- imputation_metrics(
    se_nm,
    assay_name = "cycloess",
    methods    = NULL,
    combo_methods = list(
      "Impseq+min"    = list(mar_method = "Impseq", mnar_method = "min"),
      "Impseqrob+min" = list(mar_method = "Impseqrob", mnar_method = "min")
    )
  )


  # ---- 8. From_data NA pattern (mimics real NA distribution) -----------------

  res_fd <- imputation_metrics(
    se_nm,
    assay_name = "cycloess",
    pattern    = "from_data",
    methods    = c("knn", "Impseqrob", "min", "QRILC")
  )


  # ---- 9. Export plots to PNG ------------------------------------------------

  output_dir <- "./results/imputation_metrics"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (plot_name in names(res)) {
    p <- res[[plot_name]]
    if (is.null(p) || !inherits(p, "gg")) next
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0("im_", plot_name, ".png")),
      plot     = p,
      width    = 12,
      height   = 8,
      dpi      = 150
    )
  }

}
