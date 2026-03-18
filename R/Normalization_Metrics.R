# =============================================================================
# Normalization Quality Metrics
# =============================================================================
#
# Functions for evaluating and comparing proteomics normalization methods:
#   - import_norm_matrices() : Load normalized TSV files into SummarizedExperiment
#   - normalization_metrics(): Orchestrator returning a list of ggplot2 plots
#   - nm_compute_metrics()   : Quantitative group-separation metrics (data.frame)
#
# Individual plot functions (12):
#   nm_plot_boxplot, nm_plot_density, nm_plot_pcv,
#   nm_plot_pmad, nm_plot_pev, nm_plot_pca, nm_plot_correlation,
#   nm_plot_mds, nm_plot_scatter, nm_plot_qq, nm_plot_metrics,
#   nm_plot_pc1_ranking, nm_plot_mds1_ranking
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
#   - vegan    : PERMANOVA R² (fallback: NA)
#   - cluster  : silhouette width (fallback: NA)
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# --- Null coalescing operator ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# --- Self-dir sourcing for Normalization.R ---
.self_dir <- if (sys.nframe() > 0) dirname(sys.frame(1)$ofile) else "R"

if (!exists(".norm_log2norm", mode = "function")) {
  .norm_source_path <- file.path(.self_dir, "Normalization.R")
  if (file.exists(.norm_source_path)) {
    source(.norm_source_path, local = FALSE)
  } else {
    warning("Normalization.R not found at '", .norm_source_path,
            "'. Auto-normalization will not be available.")
  }
}

# --- Benchmark methods (14, excludes "log2" which is the baseline) ---
.NM_BENCH_METHODS <- c(
  "log2Norm", "GlobalMedian", "GlobalMean", "eqmedians",
  "vsn", "max", "medianNorm", "meanNorm",
  "quantile", "Rlr", "MAD", "cycloess",
  "center_quantile", "quantile.robust"
)

# Methods that require x_raw (Grupo A) — rest use x_log2 (Grupo B)
.NM_RAW_METHODS <- c(
  "log2Norm", "GlobalMedian", "GlobalMean", "eqmedians",
  "vsn", "max", "medianNorm", "meanNorm"
)

<<<<<<< HEAD
# --- Metric direction registry (higher/lower = better) ---
# Curated set of 7 non-redundant metrics:
#   Removed: Condition_Number (no discrimination), Hopkins (misleading on log2),
#   MDS_GOF + MDS_CophCor (penalize good methods), CumVar_PC2 (redundant w/ PC1),
#   PCV_median + PEV_median (redundant w/ PMAD_median)
.NM_METRIC_DIRECTIONS <- c(
  PC1_VarPct       = "higher",
  PC1_F_ratio      = "higher",
  PERMANOVA_R2     = "higher",
  Silhouette_mean  = "higher",
  Spectral_Entropy = "lower",
  MDS1_VarPct      = "lower",
  PMAD_median      = "lower"
)

# --- Default weights: group-separation metrics dominate (75%) ---
# PC1_VarPct is the best single predictor of DE performance.
# Group-separation tier (weight 2-3) vs data-quality tier (weight 1).
.NM_DEFAULT_WEIGHTS <- c(
  PC1_VarPct       = 3,
  PC1_F_ratio      = 2,
  PERMANOVA_R2     = 2,
  Silhouette_mean  = 2,
  Spectral_Entropy = 1,
  MDS1_VarPct      = 1,
  PMAD_median      = 1
)

=======
>>>>>>> parent of 6c86c1a (Added composite_rank multiple)
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

#' PC1 variance percentage
#'
#' Percentage of total variance explained by the first principal component.
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @return Numeric scalar (0-100).
#' @keywords internal
.nm_pc1_var_pct <- function(mat) {
  if (ncol(mat) < 2 || nrow(mat) < 2) return(NA_real_)
  pca <- prcomp(t(mat), scale. = FALSE, center = TRUE)
  100 * pca$sdev[1]^2 / sum(pca$sdev^2)
}

#' PC1 F-ratio (between-group / within-group variance on PC1 scores)
#'
#' One-way ANOVA F-statistic on PC1 scores grouped by condition.
#' Higher values indicate better group separation along PC1.
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @param groups Factor or character vector of group labels.
#' @return Numeric scalar (0-Inf).
#' @keywords internal
.nm_pc1_f_ratio <- function(mat, groups) {
  if (ncol(mat) < 2 || nrow(mat) < 2) return(NA_real_)
  groups <- as.factor(groups)
  if (nlevels(groups) < 2) return(NA_real_)
  pca    <- prcomp(t(mat), scale. = FALSE, center = TRUE)
  scores <- pca$x[, 1]
  grand  <- mean(scores)
  k      <- nlevels(groups)
  n      <- length(scores)
  ss_b   <- sum(tapply(scores, groups, function(x) length(x) * (mean(x) - grand)^2))
  ss_w   <- sum(tapply(scores, groups, function(x) sum((x - mean(x))^2)))
  df_b   <- k - 1
  df_w   <- n - k
  if (df_w < 1 || ss_w == 0) return(NA_real_)
  (ss_b / df_b) / (ss_w / df_w)
}

#' PERMANOVA R² via vegan::adonis2
#'
#' Proportion of variance in Euclidean distances explained by the grouping.
#' Requires the vegan package (optional).
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @param groups Factor or character vector of group labels.
#' @return Named list with `R2` and `p_value`, or NA if vegan unavailable.
#' @keywords internal
.nm_permanova_r2 <- function(mat, groups) {
  na_result <- list(R2 = NA_real_, p_value = NA_real_)
  if (!requireNamespace("vegan", quietly = TRUE)) return(na_result)
  if (ncol(mat) < 2 || nrow(mat) < 2) return(na_result)
  groups <- as.factor(groups)
  if (nlevels(groups) < 2) return(na_result)
  d  <- dist(t(mat))
  df <- data.frame(Condition = groups)
  res <- vegan::adonis2(d ~ Condition, data = df, permutations = 999)
  list(R2 = res[["R2"]][1], p_value = res[["Pr(>F)"]][1])
}

#' Average silhouette width
#'
#' Mean silhouette width when samples are clustered by condition labels.
#' Requires the cluster package (optional).
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @param groups Factor or character vector of group labels.
#' @return Numeric scalar (-1 to 1), or NA if cluster unavailable.
#' @keywords internal
.nm_silhouette_avg <- function(mat, groups) {
  if (!requireNamespace("cluster", quietly = TRUE)) return(NA_real_)
  if (ncol(mat) < 2 || nrow(mat) < 2) return(NA_real_)
  groups <- as.factor(groups)
  if (nlevels(groups) < 2) return(NA_real_)
  d   <- dist(t(mat))
  sil <- cluster::silhouette(as.integer(groups), d)
  mean(sil[, "sil_width"])
}

#' MDS goodness-of-fit
#'
#' GOF[1] from `cmdscale()` with `eig = TRUE`: proportion of variance
#' retained in the 2D MDS projection. Uses scaled data (consistent with
#' `nm_plot_mds`).
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @return Numeric scalar (0-1).
#' @keywords internal
.nm_mds_gof <- function(mat) {
  if (ncol(mat) < 3 || nrow(mat) < 2) return(NA_real_)
  d   <- dist(scale(t(mat)))
  mds <- cmdscale(d, k = 2, eig = TRUE)
  mds$GOF[1]
}

#' Percentage of variance explained by MDS dimension 1
#'
#' Computes the percentage of variance captured by the first MDS dimension,
#' using only positive eigenvalues from classical MDS (`cmdscale(eig = TRUE)`).
#' Uses scaled data (consistent with `nm_plot_mds`).
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @return Numeric scalar (0-100).
#' @keywords internal
.nm_mds1_var_pct <- function(mat) {
  if (ncol(mat) < 3 || nrow(mat) < 2) return(NA_real_)
  d   <- dist(scale(t(mat)))
  mds <- cmdscale(d, k = 2, eig = TRUE)
  eig <- mds$eig
  pos <- eig[eig > 0]
  if (length(pos) == 0) return(NA_real_)
  100 * pos[1] / sum(pos)
}

#' MDS cophenetic correlation
#'
#' Pearson correlation between the original Euclidean distances and the
#' distances in the 2D MDS projection. Uses scaled data (consistent with
#' `nm_plot_mds`).
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @return Numeric scalar (-1 to 1).
#' @keywords internal
.nm_cophenetic_cor <- function(mat) {
  if (ncol(mat) < 3 || nrow(mat) < 2) return(NA_real_)
  d_orig <- dist(scale(t(mat)))
  mds    <- cmdscale(d_orig, k = 2)
  d_mds  <- dist(mds)
  cor(as.numeric(d_orig), as.numeric(d_mds))
}

#' Spectral entropy of PCA eigenvalues
#'
#' Normalized Shannon entropy of the eigenvalue distribution from PCA.
#' Values near 0 indicate variance concentrated in few components; values
#' near 1 indicate uniform spread across all components.
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @return Numeric scalar (0-1).
#' @keywords internal
.nm_spectral_entropy <- function(mat) {
  if (ncol(mat) < 2 || nrow(mat) < 2) return(NA_real_)
  pca <- tryCatch(prcomp(t(mat), center = TRUE, scale. = FALSE),
                  error = function(e) NULL)
  if (is.null(pca)) return(NA_real_)
  eigvals <- pca$sdev^2
  eigvals <- eigvals[eigvals > 0]
  if (length(eigvals) < 2) return(NA_real_)
  p <- eigvals / sum(eigvals)
  -sum(p * log(p)) / log(length(p))
}

#' Cumulative variance explained by PC1 and PC2
#'
#' Percentage of total variance captured by the first two principal components.
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @return Numeric scalar (0-100).
#' @keywords internal
.nm_cumvar_pc2 <- function(mat) {
  if (ncol(mat) < 2 || nrow(mat) < 2) return(NA_real_)
  pca <- tryCatch(prcomp(t(mat), center = TRUE, scale. = FALSE),
                  error = function(e) NULL)
  if (is.null(pca)) return(NA_real_)
  vars <- pca$sdev^2
  k <- min(2, length(vars))
  100 * sum(vars[seq_len(k)]) / sum(vars)
}

#' Hopkins statistic for clustering tendency
#'
#' Manual implementation of the Hopkins statistic. Values > 0.5 suggest
#' non-random clustering structure; values near 0.5 indicate uniform
#' distribution.
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @param n_sample Integer. Number of points to sample (default 10, capped at
#'   ncol - 1).
#' @return Numeric scalar (0-1).
#' @keywords internal
.nm_hopkins <- function(mat, n_sample = 10L) {
  x <- t(mat)                          # samples as rows
  n <- nrow(x)
  if (n < 3 || ncol(x) < 1) return(NA_real_)
  m <- min(n_sample, n - 1L)

  # Random reference points within data bounding box
  mins <- apply(x, 2, min)
  maxs <- apply(x, 2, max)
  rand_pts <- mapply(function(lo, hi) stats::runif(m, lo, hi),
                     mins, maxs, SIMPLIFY = TRUE)
  if (is.null(dim(rand_pts))) rand_pts <- matrix(rand_pts, nrow = m)

  # Sample m real data points (without replacement)
  set.seed(42L)
  idx <- sample.int(n, m)
  real_pts <- x[idx, , drop = FALSE]

  # Nearest-neighbour distance for random points to real data
  u <- vapply(seq_len(m), function(j) {
    min(sqrt(rowSums((sweep(x, 2, rand_pts[j, ]))^2)))
  }, numeric(1))

  # Nearest-neighbour distance for sampled real points to remaining data
  w <- vapply(seq_len(m), function(j) {
    others <- x[-idx[j], , drop = FALSE]
    min(sqrt(rowSums((sweep(others, 2, real_pts[j, ]))^2)))
  }, numeric(1))

  sum(u) / (sum(u) + sum(w))
}

#' Condition number of the covariance matrix (via PCA eigenvalues)
#'
#' Ratio of the largest to smallest non-zero eigenvalue from PCA of the
#' sample covariance matrix. High values indicate numerical instability or
#' multicollinearity.
#'
#' @param mat Numeric matrix (proteins x samples), no NAs.
#' @return Numeric scalar >= 1 (Inf if smallest eigenvalue is zero).
#' @keywords internal
.nm_condition_number <- function(mat) {
  if (ncol(mat) < 2 || nrow(mat) < 2) return(NA_real_)
  pca <- tryCatch(prcomp(t(mat), center = TRUE, scale. = FALSE),
                  error = function(e) NULL)
  if (is.null(pca)) return(NA_real_)
  eigvals <- pca$sdev^2
  eigvals <- eigvals[eigvals > 0]
  if (length(eigvals) < 2) return(NA_real_)
  max(eigvals) / min(eigvals)
}

# =============================================================================
# SECTION 1b: AUTO-NORMALIZATION DISPATCH
# =============================================================================

#' Dispatch a single normalization method on a log2 matrix
#'
#' Internal helper that calls the appropriate `.norm_*()` function.
#' For Grupo A methods, converts log2 back to raw (2^x) before calling.
#'
#' @param x_log2 Numeric matrix in log2 scale (proteins x samples)
#' @param method Character scalar: normalization method name
#' @param method_args Named list of per-method arguments
#' @return Numeric matrix in log2 scale
#' @keywords internal
.nm_dispatch_normalization <- function(x_log2, method, method_args = list()) {
  # Grupo A methods need raw (linear) scale input
  if (method %in% .NM_RAW_METHODS) {
    x_input <- 2^x_log2
  } else {
    x_input <- x_log2
  }

  x_norm <- switch(method,
    "log2Norm"        = .norm_log2norm(x_input),
    "GlobalMedian"    = .norm_ginorm(x_input),
    "GlobalMean"      = .norm_globalmean(x_input),
    "eqmedians"       = .norm_eqmedians(x_input),
    "vsn"             = .norm_vsn(x_input),
    "max"             = .norm_max(x_input),
    "medianNorm"      = .norm_mediannorm(x_input),
    "meanNorm"        = .norm_meannorm(x_input),
    "quantile"        = .norm_quantile(x_input),
    "Rlr"             = .norm_rlr(x_input),
    "MAD"             = .norm_mad(x_input),
    "cycloess"        = {
      args <- method_args[["cycloess"]] %||% list()
      limma::normalizeCyclicLoess(
        x_input,
        method     = args[["method"]]     %||% "fast",
        iterations = args[["iterations"]] %||% 3,
        span       = args[["span"]]       %||% 0.7
      )
    },
    "center_quantile" = {
      q <- (method_args[["center_quantile"]] %||% list())[["q"]] %||% 0.15
      .norm_center_quantile(x_input, q = q)
    },
    "quantile.robust" = .norm_quantile_robust(x_input),
    stop("Unknown normalization method: '", method, "'")
  )

  rownames(x_norm) <- rownames(x_log2)
  colnames(x_norm) <- colnames(x_log2)
  x_norm
}

#' Prepare a SummarizedExperiment from a spectronaut_data object
#'
#' Convenience wrapper that extracts metadata and protein data from a
#' `spectronaut_data` object (output of `preprocess_spectronaut()`), performs
#' zero-to-NA conversion, protein filtering, and returns a SE with assays
#' `"raw"` and `"log2"` — ready for `nm_run_normalizations()` or
#' `normalization_metrics(..., methods = "all")`.
#'
#' Internally calls `normalize_proteomics()` with `norm_method = "log2"`
#' (no additional normalization).
#'
#' @param preprocessing `spectronaut_data` list from `preprocess_spectronaut()`.
#' @param min_reps Minimum replicates with non-NA values per group for protein
#'   filtering. If NULL, auto-computed as half the smallest group. Default `NULL`.
#' @param min_groups Minimum groups meeting `min_reps` (default: 1).
#' @param covariate_df Optional covariate data.frame for paired designs
#'   (must contain a `Column` column). Default `NULL`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @return SummarizedExperiment with assays `"raw"` and `"log2"`.
#'
#' @examples
#' \dontrun{
#' source("R/Normalization_Metrics.R")
#' se <- nm_prepare_se(preprocessing, min_reps = 3)
#' plots <- normalization_metrics(se, methods = "all")
#' }
#' @export
nm_prepare_se <- function(preprocessing,
                          min_reps     = NULL,
                          min_groups   = 1,
                          covariate_df = NULL,
                          verbose      = TRUE) {

  if (!inherits(preprocessing, "spectronaut_data"))
    stop("'preprocessing' must be a spectronaut_data object ",
         "(output of preprocess_spectronaut()).")

  # Source Processing.R for .prepare_metadata / .prepare_protein_data
  if (!exists(".prepare_metadata", mode = "function")) {
    proc_path <- file.path(.self_dir, "Processing.R")
    if (file.exists(proc_path)) {
      source(proc_path, local = FALSE)
    } else {
      stop("Processing.R not found at '", proc_path,
           "'. Required for nm_prepare_se().")
    }
  }

  metadata     <- .prepare_metadata(preprocessing, covariate_df = covariate_df)
  protein_data <- .prepare_protein_data(preprocessing)

  norm_result <- normalize_proteomics(
    data       = protein_data,
    metadata   = metadata,
    min_reps   = min_reps,
    min_groups = min_groups,
    norm_method = "log2",
    verbose     = verbose
  )

  norm_result$se
}

#' Run multiple normalization methods from a baseline assay
#'
#' Takes a SummarizedExperiment with a log2-scale assay and applies each
#' requested normalization method, returning a new SE with one assay per method.
#'
#' @param se SummarizedExperiment with at least one log2-scale assay.
#' @param assay_name Name of the baseline assay to normalize from.
#'   Default `"log2"`.
#' @param methods Character vector of method names, or `"all"` for all 14
#'   benchmark methods. Default `"all"`.
#' @param method_args Named list of per-method arguments. E.g.
#'   `list(cycloess = list(method = "fast", span = 0.8))`.
#' @param include_baseline Logical. Include the baseline assay in the output SE.
#'   Default `TRUE`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @return SummarizedExperiment with one assay per successfully normalized method
#'   (plus baseline if `include_baseline = TRUE`).
#'
#' @examples
#' \dontrun{
#' se_bench <- nm_run_normalizations(se, assay_name = "log2", methods = "all")
#' SummarizedExperiment::assayNames(se_bench)
#' }
#' @export
nm_run_normalizations <- function(se,
                                  assay_name       = "log2",
                                  methods          = "all",
                                  method_args      = list(),
                                  include_baseline = TRUE,
                                  verbose          = TRUE) {

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("Package 'SummarizedExperiment' is required.")

  # Validate baseline assay exists
  all_assays <- SummarizedExperiment::assayNames(se)
  if (!assay_name %in% all_assays)
    stop("Assay '", assay_name, "' not found in SE. Available: ",
         paste(all_assays, collapse = ", "))

  # Resolve methods
  if (identical(methods, "all")) {
    methods <- .NM_BENCH_METHODS
  } else {
    unknown <- setdiff(methods, .NM_BENCH_METHODS)
    if (length(unknown) > 0)
      warning("Unknown method(s) ignored: ", paste(unknown, collapse = ", "))
    methods <- intersect(methods, .NM_BENCH_METHODS)
    if (length(methods) == 0)
      stop("No valid normalization methods provided.")
  }

  x_log2 <- SummarizedExperiment::assay(se, assay_name)

  # log2Norm is identical to log2 baseline — drop it when baseline is included
  if (include_baseline && "log2Norm" %in% methods) {
    methods <- setdiff(methods, "log2Norm")
    if (verbose) message("  Skipping 'log2Norm' (identical to baseline '",
                         assay_name, "')")
  }

  # Build assay list
  assay_list <- list()
  if (include_baseline) assay_list[[assay_name]] <- x_log2

  n_ok   <- 0L
  n_fail <- 0L

  for (m in methods) {
    if (verbose) message("  Normalizing: ", m, " ...")
    result <- tryCatch(
      .nm_dispatch_normalization(x_log2, m, method_args),
      error = function(e) {
        warning("Method '", m, "' failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(result)) {
      assay_list[[m]] <- result
      n_ok <- n_ok + 1L
    } else {
      n_fail <- n_fail + 1L
    }
  }

  if (verbose)
    message("nm_run_normalizations: ", n_ok, " method(s) OK",
            if (n_fail > 0) paste0(", ", n_fail, " failed") else ".")

  # Build new SE
  se_out <- SummarizedExperiment::SummarizedExperiment(
    assays  = assay_list,
    colData = SummarizedExperiment::colData(se),
    rowData = if (nrow(SummarizedExperiment::rowData(se)) > 0)
                SummarizedExperiment::rowData(se) else NULL
  )
  se_out
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
#' @param pca_scales Facet scaling: `"free"` (default) allows independent axes
#'   per method; `"fixed"` uses shared axes to compare separation magnitude.
#' @return ggplot object.
#' @export
nm_plot_pca <- function(se, assay_names = NULL,
                        condition_col = "Condition",
                        pca_scales = c("free", "fixed"), ...) {
  pca_scales <- match.arg(pca_scales)
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
    ggplot2::facet_wrap(~ Facet, ncol = 2, scales = pca_scales) +
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
#' @param mds_scales Facet scaling: `"free"` (default) allows independent axes
#'   per panel; `"fixed"` forces shared axes for easier cross-method comparison.
#' @return ggplot object.
#' @export
nm_plot_mds <- function(se, assay_names = NULL,
                        condition_col = "Condition",
                        mds_scales = c("free", "fixed"), ...) {
  mds_scales <- match.arg(mds_scales)
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
    ggplot2::facet_wrap(~ Method, ncol = 2, scales = mds_scales) +
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

# --------------------------------------------------------------------------
# 11. Metrics — quantitative group-separation metrics
# --------------------------------------------------------------------------

#' Compute quantitative group-separation metrics per normalization method
#'
#' Iterates over assays in a SummarizedExperiment and computes ten metrics
#' that quantify how well the normalization separates sample groups.
#'
#' @inheritParams nm_plot_boxplot
#' @return A `data.frame` with one row per method and columns:
#'   `Method`, `PC1_VarPct`, `PC1_F_ratio`, `PERMANOVA_R2`, `PERMANOVA_pval`,
#'   `Silhouette_mean`, `MDS_GOF`, `MDS_CophCor`, `Spectral_Entropy`,
#'   `CumVar_PC2`, `Hopkins`, `Condition_Number`.
#'
#' @details
#' - **PC1_VarPct**: % variance explained by PC1 (higher = more structure).
#' - **PC1_F_ratio**: ANOVA F-ratio on PC1 scores (higher = better separation).
#' - **PERMANOVA_R2**: Proportion of variance explained by grouping (requires
#'   vegan; NA if not installed).
#' - **Silhouette_mean**: Mean silhouette width (requires cluster; NA if not
#'   installed).
#' - **MDS_GOF**: Goodness-of-fit of 2D MDS projection. Note: in benchmarking
#'   context, *lower* values tend to indicate better normalization (a good method
#'   preserves multi-dimensional biological variability that 2D cannot capture).
#' - **MDS_CophCor**: Cophenetic correlation between original and MDS distances.
#'   Same caveat as MDS_GOF: *lower* values tend to correlate with better
#'   benchmarking performance.
#' - **Spectral_Entropy**: Normalized Shannon entropy of PCA eigenvalues (0-1).
#'   Values near 0 = variance concentrated in few PCs; near 1 = uniform spread.
#' - **CumVar_PC2**: Cumulative % variance explained by PC1 + PC2 (0-100).
#' - **Hopkins**: Hopkins statistic for clustering tendency (0-1). Values > 0.5
#'   suggest non-random cluster structure.
#' - **Condition_Number**: Ratio of largest to smallest PCA eigenvalue (>= 1).
#'   High values indicate multicollinearity or numerical instability.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' metrics_df <- nm_compute_metrics(se_nm)
#' print(metrics_df)
#' }
#' @export
nm_compute_metrics <- function(se, assay_names = NULL,
                               condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)

  # Inform about optional packages
  if (!requireNamespace("vegan", quietly = TRUE))
    message("nm_compute_metrics: 'vegan' not installed — PERMANOVA columns will be NA.")
  if (!requireNamespace("cluster", quietly = TRUE))
    message("nm_compute_metrics: 'cluster' not installed — Silhouette column will be NA.")

  rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat    <- SummarizedExperiment::assay(se, assay_names[i])
    mat_ok <- mat[complete.cases(mat), ]
    groups <- condition

    perm <- .nm_permanova_r2(mat_ok, groups)

    rows[[i]] <- data.frame(
      Method          = assay_names[i],
      PC1_VarPct      = .nm_pc1_var_pct(mat_ok),
      PC1_F_ratio     = .nm_pc1_f_ratio(mat_ok, groups),
      PERMANOVA_R2    = perm$R2,
      PERMANOVA_pval  = perm$p_value,
      Silhouette_mean = .nm_silhouette_avg(mat_ok, groups),
      MDS_GOF         = .nm_mds_gof(mat_ok),
      MDS_CophCor        = .nm_cophenetic_cor(mat_ok),
      Spectral_Entropy   = .nm_spectral_entropy(mat_ok),
      CumVar_PC2         = .nm_cumvar_pc2(mat_ok),
      Hopkins            = .nm_hopkins(mat_ok),
      Condition_Number   = .nm_condition_number(mat_ok),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Bar chart of group-separation metrics per normalization method
#'
#' Calls `nm_compute_metrics()` internally and produces a faceted bar chart
#' (one facet per metric, free y-scales) with labeled values.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' nm_plot_metrics(se_nm)
#' }
#' @export
nm_plot_metrics <- function(se, assay_names = NULL,
                            condition_col = "Condition", ...) {
  assay_names <- .nm_assay_names(se, assay_names)
  metrics_df  <- nm_compute_metrics(se, assay_names, condition_col)
  col_vector  <- .nm_prone_colors(length(assay_names))

  # Pivot to long format (exclude PERMANOVA_pval from plot)
  value_cols <- c("PC1_VarPct", "PC1_F_ratio", "PERMANOVA_R2",
                  "Silhouette_mean", "MDS_GOF", "MDS_CophCor",
                  "Spectral_Entropy", "CumVar_PC2", "Hopkins",
                  "Condition_Number")
  long_df <- tidyr::pivot_longer(
    metrics_df[, c("Method", value_cols)],
    cols      = tidyr::all_of(value_cols),
    names_to  = "Metric",
    values_to = "Value"
  )
  long_df$Method <- factor(long_df$Method, levels = assay_names)
  long_df$Metric <- factor(long_df$Metric, levels = value_cols)

  ggplot2::ggplot(long_df,
    ggplot2::aes(x = Method, y = Value, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_label(
      ggplot2::aes(label = ifelse(is.na(Value), "NA",
                                  sprintf("%.2f", Value))),
      size = 2.5, fill = "white", linewidth = 0.2, na.rm = TRUE) +
    ggplot2::facet_wrap(~ Metric, ncol = 2, scales = "free_y") +
    ggplot2::scale_fill_manual(values = col_vector) +
    ggplot2::labs(
      title = "Group-Separation Metrics per Normalization Method",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text    = ggplot2::element_text(face = "bold")
    )
}

# --------------------------------------------------------------------------
# 12. PC1 Variance Ranking
# --------------------------------------------------------------------------

#' Rank normalization methods by PC1 variance explained
#'
#' Computes the percentage of total variance captured by PC1 for each assay
#' in a SummarizedExperiment and returns a data.frame sorted in descending
#' order. Higher PC1 variance generally indicates stronger group separation.
#'
#' @inheritParams nm_plot_boxplot
#' @return A `data.frame` with columns `Method`, `PC1_VarPct`, `Rank`,
#'   ordered by `PC1_VarPct` descending.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' nm_rank_pc1(se_nm)
#' }
#' @export
nm_rank_pc1 <- function(se, assay_names = NULL, condition_col = "Condition",
                        verbose = TRUE) {
  assay_names <- .nm_assay_names(se, assay_names)

  pct <- vapply(assay_names, function(nm) {
    mat    <- SummarizedExperiment::assay(se, nm)
    mat_ok <- mat[complete.cases(mat), ]
    .nm_pc1_var_pct(mat_ok)
  }, numeric(1))

  df <- data.frame(
    Method     = assay_names,
    PC1_VarPct = pct,
    stringsAsFactors = FALSE
  )
  df <- df[order(-df$PC1_VarPct), ]
  df$Rank <- seq_len(nrow(df))
  rownames(df) <- NULL

  df
}

#' Horizontal bar chart of PC1 variance ranking
#'
#' Produces a horizontal bar chart with normalization methods ordered by
#' descending PC1 variance percentage. Uses the same PRONE-style palette
#' as other `nm_plot_*()` functions.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' nm_plot_pc1_ranking(se_nm)
#' }
#' @export
nm_plot_pc1_ranking <- function(se, assay_names = NULL,
                                condition_col = "Condition", ...) {
  rank_df    <- nm_rank_pc1(se, assay_names, condition_col, verbose = FALSE)
  col_vector <- .nm_prone_colors(nrow(rank_df))

  # Order factor by PC1_VarPct descending (bottom-to-top in coord_flip)
  rank_df$Method <- factor(rank_df$Method,
                           levels = rev(rank_df$Method))

  ggplot2::ggplot(rank_df,
    ggplot2::aes(x = Method, y = PC1_VarPct, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", PC1_VarPct)),
      hjust = -0.1, size = 3.2) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = rep_len(col_vector, nrow(rank_df))) +
    ggplot2::labs(
      title = "PC1 Variance Explained per Normalization Method",
      x     = NULL,
      y     = "PC1 Variance (%)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::expand_limits(y = max(rank_df$PC1_VarPct, na.rm = TRUE) * 1.08)
}

#' Rank normalization methods by MDS1 variance explained
#'
#' For each assay in `se`, computes the percentage of variance captured by the
#' first MDS dimension (using only positive eigenvalues from classical MDS) and
#' returns a ranking ordered by ascending MDS1 variance (lower = better).
#'
#' @inheritParams nm_plot_boxplot
#' @return A `data.frame` with columns `Method`, `MDS1_VarPct`, `Rank`,
#'   ordered by `MDS1_VarPct` descending.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' nm_rank_mds1(se_nm)
#' }
#' @export
nm_rank_mds1 <- function(se, assay_names = NULL, condition_col = "Condition",
                         verbose = TRUE) {
  assay_names <- .nm_assay_names(se, assay_names)

  pct <- vapply(assay_names, function(nm) {
    mat    <- SummarizedExperiment::assay(se, nm)
    mat_ok <- mat[complete.cases(mat), ]
    .nm_mds1_var_pct(mat_ok)
  }, numeric(1))

  df <- data.frame(
    Method      = assay_names,
    MDS1_VarPct = pct,
    stringsAsFactors = FALSE
  )
  df <- df[order(df$MDS1_VarPct), ]
  df$Rank <- seq_len(nrow(df))
  rownames(df) <- NULL

  df
}

#' Horizontal bar chart of MDS1 variance ranking
#'
#' Produces a horizontal bar chart with normalization methods ordered by
#' ascending MDS1 variance percentage (lower = better normalization).
#' Uses the same PRONE-style palette as other `nm_plot_*()` functions.
#'
#' @inheritParams nm_plot_boxplot
#' @return ggplot object.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' nm_plot_mds1_ranking(se_nm)
#' }
#' @export
nm_plot_mds1_ranking <- function(se, assay_names = NULL,
                                 condition_col = "Condition", ...) {
  rank_df    <- nm_rank_mds1(se, assay_names, condition_col, verbose = FALSE)
  col_vector <- .nm_prone_colors(nrow(rank_df))

  # Order factor by MDS1_VarPct ascending (best at top in coord_flip)
  rank_df$Method <- factor(rank_df$Method,
                           levels = rev(rank_df$Method))

  ggplot2::ggplot(rank_df,
    ggplot2::aes(x = Method, y = MDS1_VarPct, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", MDS1_VarPct)),
      hjust = -0.1, size = 3.2) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = rep_len(col_vector, nrow(rank_df))) +
    ggplot2::labs(
      title = "MDS1 Variance Explained per Normalization Method",
      x     = NULL,
      y     = "MDS1 Variance (%)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::expand_limits(y = max(rank_df$MDS1_VarPct, na.rm = TRUE) * 1.08)
}

<<<<<<< HEAD
# --------------------------------------------------------------------------
# 14. Composite Ranking (rank-aggregation across all metrics)
# --------------------------------------------------------------------------

#' Composite ranking of normalization methods across 7 curated quality metrics
#'
#' Combines 7 non-redundant metrics (5 from `nm_compute_metrics()` plus
#' MDS1_VarPct and PMAD_median) into a single rank-aggregation table. For each
#' metric, methods are ranked according to `.NM_METRIC_DIRECTIONS`
#' (higher-is-better or lower-is-better). Each metric is min-max normalized
#' to [0, 1] (1 = best) and `Score_Mean` is the weighted mean of scores.
#' Per-metric ranks and `Rank_Mean` are also included for reference.
#'
#' @inheritParams nm_plot_boxplot
#' @param weights Named numeric vector of metric weights. Names must match
#'   metric column names. `"default"` (the default) uses `.NM_DEFAULT_WEIGHTS`
#'   which prioritizes group-separation metrics (PC1_VarPct = 3, F_ratio /
#'   PERMANOVA / Silhouette = 2, rest = 1). Pass `NULL` for equal weights, or
#'   a custom named vector to override.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @return A `data.frame` with columns: `Method`, 7 value columns,
#'   7 `*_Rank` columns, `Rank_Mean`, 7 `*_Score` columns (min-max
#'   normalized 0-1 where 1 = best), `Score_Mean` (weighted), `Composite_Rank`
#'   (ordered by `Score_Mean` descending).
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' nm_rank_composite(se_nm)
#' }
#' @export
nm_rank_composite <- function(se, assay_names = NULL,
                              condition_col = "Condition",
                              weights = "default",
                              verbose = TRUE) {
  assay_names <- .nm_assay_names(se, assay_names)
  condition   <- .nm_condition(se, condition_col)

  # Resolve weights: "default" → .NM_DEFAULT_WEIGHTS, NULL → equal weights
  if (is.character(weights) && identical(weights, "default")) {
    weights <- .NM_DEFAULT_WEIGHTS
  }

  # --- Step 1: Base metrics from nm_compute_metrics() (keep 5 of 11) ---
  base_df <- nm_compute_metrics(se, assay_names, condition_col)
  base_keep <- c("Method", "PC1_VarPct", "PC1_F_ratio", "PERMANOVA_R2",
                 "Silhouette_mean", "Spectral_Entropy")
  base_df <- base_df[, intersect(base_keep, colnames(base_df)), drop = FALSE]

  # --- Step 2: Additional metrics (MDS1_VarPct, PMAD_median) ---
  extra_rows <- vector("list", length(assay_names))
  for (i in seq_along(assay_names)) {
    mat    <- SummarizedExperiment::assay(se, assay_names[i])
    mat_ok <- mat[complete.cases(mat), ]
    groups <- condition

    extra_rows[[i]] <- data.frame(
      Method      = assay_names[i],
      MDS1_VarPct = .nm_mds1_var_pct(mat_ok),
      PMAD_median = median(.nm_pmad(mat_ok, groups), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  extra_df <- do.call(rbind, extra_rows)

  # --- Step 3: Merge ---
  full_df <- merge(base_df, extra_df, by = "Method", sort = FALSE)

  # --- Step 4: Determine metric columns ---
  metric_cols <- setdiff(colnames(full_df), "Method")

  # Drop all-NA columns (optional deps missing, e.g. vegan, cluster)
  all_na <- vapply(metric_cols, function(m) all(is.na(full_df[[m]])), logical(1))
  if (any(all_na)) {
    if (verbose) message("nm_rank_composite: dropping all-NA metric(s): ",
                         paste(metric_cols[all_na], collapse = ", "))
    metric_cols <- metric_cols[!all_na]
  }

  if (length(metric_cols) == 0) stop("No metrics available for ranking.")

  # --- Step 5: Compute ranks and scores per metric ---
  directions <- .NM_METRIC_DIRECTIONS
  rank_df  <- full_df[, "Method", drop = FALSE]
  score_df <- full_df[, "Method", drop = FALSE]

  for (m in metric_cols) {
    x   <- full_df[[m]]
    dir <- directions[m]
    if (is.na(dir) || is.null(dir)) dir <- "higher"

    # Ranks (kept for reference)
    if (dir == "higher") {
      rank_df[[paste0(m, "_Rank")]] <- rank(-x, ties.method = "average",
                                             na.last = "keep")
    } else {
      rank_df[[paste0(m, "_Rank")]] <- rank(x, ties.method = "average",
                                             na.last = "keep")
    }

    # Min-max score [0, 1] where 1 = best
    rng <- range(x, na.rm = TRUE)
    if (rng[2] - rng[1] < .Machine$double.eps) {
      score_df[[paste0(m, "_Score")]] <- ifelse(is.na(x), NA_real_, 0.5)
    } else if (dir == "higher") {
      score_df[[paste0(m, "_Score")]] <- (x - rng[1]) / (rng[2] - rng[1])
    } else {
      score_df[[paste0(m, "_Score")]] <- (rng[2] - x) / (rng[2] - rng[1])
    }
  }

  # --- Step 6: Weighted Score_Mean ---
  score_cols <- paste0(metric_cols, "_Score")
  score_mat  <- as.matrix(score_df[, score_cols, drop = FALSE])

  if (!is.null(weights)) {
    w <- weights[metric_cols]
    w[is.na(w)] <- 1
  } else {
    w <- rep(1, length(metric_cols))
    names(w) <- metric_cols
  }

  score_df$Score_Mean <- apply(score_mat, 1, function(s) {
    ok <- !is.na(s)
    if (!any(ok)) return(NA_real_)
    stats::weighted.mean(s[ok], w[ok])
  })

  # Also compute Rank_Mean for reference
  rank_cols <- paste0(metric_cols, "_Rank")
  rank_mat  <- as.matrix(rank_df[, rank_cols, drop = FALSE])
  rank_df$Rank_Mean <- apply(rank_mat, 1, function(r) {
    ok <- !is.na(r)
    if (!any(ok)) return(NA_real_)
    stats::weighted.mean(r[ok], w[ok])
  })

  # --- Step 7: Combine values + ranks + scores, sort by Score_Mean ---
  out <- merge(full_df[, c("Method", metric_cols)], rank_df, by = "Method",
               sort = FALSE)
  out <- merge(out, score_df[, c("Method", score_cols, "Score_Mean")],
               by = "Method", sort = FALSE)
  out$Rank_Mean <- rank_df$Rank_Mean[match(out$Method, rank_df$Method)]
  out <- out[order(-out$Score_Mean), ]
  out$Composite_Rank <- seq_len(nrow(out))
  rownames(out) <- NULL

  if (verbose) {
    message("nm_rank_composite: Rank 1 = ", out$Method[1],
            " (Score_Mean = ", sprintf("%.4f", out$Score_Mean[1]),
            ", Rank_Mean = ", sprintf("%.2f", out$Rank_Mean[1]), ")")
  }

  out
}

#' Horizontal bar chart of composite normalization ranking
#'
#' Produces a horizontal bar chart with normalization methods ordered by
#' descending `Score_Mean` (higher = better). Analogous to
#' `nm_plot_pc1_ranking()` and `nm_plot_mds1_ranking()`.
#'
#' @inheritParams nm_rank_composite
#' @param ... Additional arguments (currently unused).
#' @return ggplot object.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' nm_plot_composite_ranking(se_nm)
#' }
#' @export
nm_plot_composite_ranking <- function(se, assay_names = NULL,
                                      condition_col = "Condition",
                                      weights = "default", ...) {
  comp_df    <- nm_rank_composite(se, assay_names, condition_col,
                                  weights = weights,
                                  verbose = FALSE)
  col_vector <- .nm_prone_colors(nrow(comp_df))

  # Order factor: best (highest Score_Mean) at top in coord_flip
  comp_df$Method <- factor(comp_df$Method,
                           levels = rev(comp_df$Method))

  ggplot2::ggplot(comp_df,
    ggplot2::aes(x = Method, y = Score_Mean, fill = Method)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", Score_Mean)),
      hjust = -0.1, size = 3.2) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = rep_len(col_vector, nrow(comp_df))) +
    ggplot2::labs(
      title    = "Composite Normalization Ranking",
      subtitle = "Higher weighted score = better overall performance",
      x = NULL,
      y = "Weighted Score Mean"
    ) +
    ggplot2::theme_bw() +
    ggplot2::expand_limits(y = max(comp_df$Score_Mean, na.rm = TRUE) * 1.08)
}

#' Heatmap of per-metric scores for normalization methods
#'
#' Tile heatmap showing the normalized score (0-1, 1 = best) each
#' normalization method achieved in each quality metric. Rows = methods
#' (ordered by `Score_Mean`, best at top), columns = metrics. Follows the
#' same pattern as `im_plot_ranking()` in `Imputation_Metrics.R`.
#'
#' @inheritParams nm_rank_composite
#' @param ... Additional arguments (currently unused).
#' @return ggplot object.
#'
#' @examples
#' \dontrun{
#' se_nm <- import_norm_matrices("./results", "./data/metadata.tsv")
#' nm_plot_composite_heatmap(se_nm)
#' }
#' @export
nm_plot_composite_heatmap <- function(se, assay_names = NULL,
                                      condition_col = "Condition",
                                      weights = "default", ...) {
  comp_df <- nm_rank_composite(se, assay_names, condition_col,
                               weights = weights,
                               verbose = FALSE)
  method_order <- comp_df$Method

  # Collect score columns + Score_Mean
  score_cols <- grep("_Score$|^Score_Mean$", colnames(comp_df), value = TRUE)
  score_df   <- comp_df[, c("Method", score_cols)]

  long_df <- tidyr::pivot_longer(
    score_df,
    cols      = tidyr::all_of(score_cols),
    names_to  = "Metric",
    values_to = "Score"
  )

  long_df$Method <- factor(long_df$Method, levels = rev(method_order))
  long_df$Metric <- factor(long_df$Metric, levels = score_cols)

  ggplot2::ggplot(long_df,
    ggplot2::aes(x = Metric, y = Method, fill = Score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(Score), "NA", sprintf("%.2f", Score))),
      size = 3, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
      midpoint = 0.5,
      na.value = "grey80",
      name = "Score") +
    ggplot2::labs(
      title    = "Composite Normalization Ranking — Per-Metric Scores",
      subtitle = "Higher score (blue) = better performance (0-1 normalized)",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y  = ggplot2::element_text(size = 9),
      plot.title   = ggplot2::element_text(face = "bold", size = 12),
      panel.grid   = ggplot2::element_blank()
    )
}

=======
>>>>>>> parent of 6c86c1a (Added composite_rank multiple)
# =============================================================================
# SECTION 4: EXPORT HELPER + MAIN ORCHESTRATOR
# =============================================================================

# Internal: export tables and plots to output_dir
.nm_export_results <- function(result, output_dir, export_plots, export_tables,
                               width, height, dpi, verbose) {
  if (is.null(output_dir)) return(invisible(NULL))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Tablas
  if (export_tables) {
    if (!is.null(result$metrics_table))
      utils::write.table(result$metrics_table,
                         file.path(output_dir, "nm_metrics_table.tsv"),
                         sep = "\t", row.names = FALSE, quote = FALSE)
    if (!is.null(result$pc1_rank))
      utils::write.table(result$pc1_rank,
                         file.path(output_dir, "nm_pc1_rank.tsv"),
                         sep = "\t", row.names = FALSE, quote = FALSE)
    if (!is.null(result$mds1_rank))
      utils::write.table(result$mds1_rank,
                         file.path(output_dir, "nm_mds1_rank.tsv"),
                         sep = "\t", row.names = FALSE, quote = FALSE)
    if (verbose) message("Exported tables to: ", output_dir)
  }

  # Plots
  if (export_plots) {
    n_saved <- 0L
    for (nm in names(result)) {
      obj <- result[[nm]]
      if (!is.null(obj) && inherits(obj, "gg")) {
        ggplot2::ggsave(file.path(output_dir, paste0("nm_", nm, ".png")),
                        plot = obj, width = width, height = height, dpi = dpi)
        n_saved <- n_saved + 1L
      }
    }
    if (verbose) message("Exported ", n_saved, " plot(s) to: ", output_dir)
  }
  invisible(NULL)
}

#' Generate quality metric plots for normalization comparison
#'
#' Calls all (or selected) `nm_plot_*()` functions on a SummarizedExperiment
#' produced by `import_norm_matrices()`. Each plot is wrapped in `tryCatch()`
#' so that a failure in one plot does not abort the entire run.
#'
#' When `methods` is provided, the function first runs
#' `nm_run_normalizations()` to automatically apply the requested normalization
#' methods from the `base_assay`, then computes metrics and generates plots on
#' the resulting multi-assay SE.
#'
#' @param se SummarizedExperiment from `import_norm_matrices()`.
#' @param assay_names Character vector of assay names to include.
#'   NULL (default) = all assays.
#' @param condition_col Column name in `colData(se)` with condition labels.
#'   Default `"Condition"`.
#' @param plots Character vector of plot names to generate, or `"all"` (default).
#'   Valid names: `"boxplot"`, `"density"`, `"pcv"`, `"pmad"`, `"pev"`,
#'   `"pca"`, `"correlation"`, `"mds"`, `"scatter"`, `"qq"`, `"metrics"`,
#'   `"pc1_ranking"`. When `pca_scales = "both"`, `"pca"` expands to
#'   `"pca_free"` + `"pca_fixed"`. When `mds_scales = "both"`, `"mds"` expands
#'   to `"mds_free"` + `"mds_fixed"`.
#' @param cor_method Correlation method for `nm_plot_correlation()`.
#'   Default `"pearson"`.
#' @param pca_scales Facet scaling for PCA plot: `"free"` (default),
#'   `"fixed"`, or `"both"` to generate and export both variants
#'   (`pca_free` and `pca_fixed` in the returned list).
#' @param mds_scales Facet scaling for MDS plot: `"free"` (default),
#'   `"fixed"`, or `"both"` to generate and export both variants
#'   (`mds_free` and `mds_fixed` in the returned list).
#' @param methods Character vector of normalization method names to
#'   auto-benchmark, `"all"` for all 14 methods, or NULL (default) to skip
#'   auto-normalization and use existing assays.
#' @param method_args Named list of per-method arguments forwarded to
#'   `nm_run_normalizations()`. E.g.
#'   `list(cycloess = list(method = "fast", span = 0.8))`.
#' @param base_assay Name of the baseline assay to normalize from when using
#'   auto-normalization. Default `"log2"`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @param output_dir Character. Path to export directory. If `NULL` (default),
#'   no files are exported. When set, tables (TSV) and/or plots (PNG) are
#'   saved to this directory (created if needed).
#' @param export_plots Logical. Export plots as PNG when `output_dir` is set.
#'   Default `TRUE`.
#' @param export_tables Logical. Export tables as TSV when `output_dir` is set.
#'   Default `TRUE`.
#' @param plot_width Numeric. Width in inches for exported plots. Default `12`.
#' @param plot_height Numeric. Height in inches for exported plots. Default `8`.
#' @param plot_dpi Numeric. Resolution for exported plots. Default `150`.
#' @return Named list of ggplot objects (or NULL for failed plots), plus
#'   `metrics_table`: a `data.frame` from `nm_compute_metrics()` and
#'   `pc1_rank`: a `data.frame` from `nm_rank_pc1()` (both always computed
#'   regardless of `plots` selection).
#'
#' @examples
#' \dontrun{
#' # --- Classic workflow (from pre-computed TSVs) ---
#' se_nm  <- import_norm_matrices("./results", "./data/metadata.tsv")
#' plots  <- normalization_metrics(se_nm)
#' plots$scatter
#' plots$pca
#'
#' # --- Auto-benchmark: all methods from a single SE ---
#' plots <- normalization_metrics(se, methods = "all", base_assay = "log2")
#'
#' # --- Auto-benchmark: specific methods ---
#' plots <- normalization_metrics(se, methods = c("cycloess", "MAD", "vsn"))
#'
#' # --- With custom parameters ---
#' plots <- normalization_metrics(se, methods = "all",
#'   method_args = list(cycloess = list(method = "fast", span = 0.8)))
#'
#' # --- With auto-export ---
#' plots <- normalization_metrics(se_nm, output_dir = "output/norm_metrics")
#' # Creates output/norm_metrics/ with nm_*.tsv and nm_*.png
#' }
#' @export
normalization_metrics <- function(se,
                                  assay_names   = NULL,
                                  condition_col = "Condition",
                                  plots         = "all",
                                  cor_method    = "pearson",
                                  pca_scales    = c("free", "fixed", "both"),
                                  mds_scales    = c("free", "fixed", "both"),
                                  methods       = NULL,
                                  method_args   = list(),
                                  base_assay    = "log2",
                                  verbose       = TRUE,
                                  output_dir    = NULL,
                                  export_plots  = TRUE,
                                  export_tables = TRUE,
                                  plot_width    = 12,
                                  plot_height   = 8,
                                  plot_dpi      = 150) {
  pca_scales <- match.arg(pca_scales)
  mds_scales <- match.arg(mds_scales)

  # --- Required packages check ---
  for (pkg in c("ggplot2", "dplyr", "tidyr", "SummarizedExperiment", "S4Vectors")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for normalization_metrics().")
  }

  # --- Auto-normalization (if methods requested) ---
  if (!is.null(methods)) {
    if (verbose) message("Auto-normalizing from assay '", base_assay, "' ...")
    se <- nm_run_normalizations(se, assay_name = base_assay,
                                methods = methods, method_args = method_args,
                                include_baseline = TRUE, verbose = verbose)
    if (is.null(assay_names))
      assay_names <- SummarizedExperiment::assayNames(se)
  }

  # --- Resolve assay names ---
  all_assays  <- SummarizedExperiment::assayNames(se)
  assay_names <- assay_names %||% all_assays
  missing_a   <- setdiff(assay_names, all_assays)
  if (length(missing_a) > 0)
    stop("Assay(s) not found in SE: ", paste(missing_a, collapse = ", "))

  # --- Plot registry ---
  all_plot_names <- c("boxplot", "density", "pcv", "pmad", "pev",
                      "pca", "correlation", "mds", "scatter", "qq",
                      "metrics", "pc1_ranking", "mds1_ranking")

  # When pca_scales == "both", expand "pca" into "pca_free" + "pca_fixed"
  if (pca_scales == "both") {
    all_plot_names <- c(setdiff(all_plot_names, "pca"), "pca_free", "pca_fixed")
  }
  # When mds_scales == "both", expand "mds" into "mds_free" + "mds_fixed"
  if (mds_scales == "both") {
    all_plot_names <- c(setdiff(all_plot_names, "mds"), "mds_free", "mds_fixed")
  }

  plot_fns <- list(
    boxplot     = function() nm_plot_boxplot(se, assay_names, condition_col),
    density     = function() nm_plot_density(se, assay_names, condition_col),
    pcv         = function() nm_plot_pcv(se, assay_names, condition_col),
    pmad        = function() nm_plot_pmad(se, assay_names, condition_col),
    pev         = function() nm_plot_pev(se, assay_names, condition_col),
    correlation = function() nm_plot_correlation(se, assay_names, condition_col,
                                                 cor_method = cor_method),
    #mds — set conditionally below
    scatter     = function() nm_plot_scatter(se, assay_names, condition_col),
    qq          = function() nm_plot_qq(se, assay_names, condition_col),
    metrics     = function() nm_plot_metrics(se, assay_names, condition_col),
    pc1_ranking  = function() nm_plot_pc1_ranking(se, assay_names, condition_col),
    mds1_ranking = function() nm_plot_mds1_ranking(se, assay_names, condition_col)
  )
  if (pca_scales == "both") {
    plot_fns$pca_free  <- function() nm_plot_pca(se, assay_names, condition_col,
                                                  pca_scales = "free")
    plot_fns$pca_fixed <- function() nm_plot_pca(se, assay_names, condition_col,
                                                  pca_scales = "fixed")
  } else {
    plot_fns$pca <- function() nm_plot_pca(se, assay_names, condition_col,
                                            pca_scales = pca_scales)
  }
  if (mds_scales == "both") {
    plot_fns$mds_free  <- function() nm_plot_mds(se, assay_names, condition_col,
                                                  mds_scales = "free")
    plot_fns$mds_fixed <- function() nm_plot_mds(se, assay_names, condition_col,
                                                  mds_scales = "fixed")
  } else {
    plot_fns$mds <- function() nm_plot_mds(se, assay_names, condition_col,
                                            mds_scales = mds_scales)
  }

  # --- Determine which plots to run ---
  if (identical(plots, "all")) {
    selected <- all_plot_names
  } else {
    # Expand "pca" → "pca_free" + "pca_fixed" when pca_scales == "both"
    if (pca_scales == "both" && "pca" %in% plots)
      plots <- c(setdiff(plots, "pca"), "pca_free", "pca_fixed")
    # Expand "mds" → "mds_free" + "mds_fixed" when mds_scales == "both"
    if (mds_scales == "both" && "mds" %in% plots)
      plots <- c(setdiff(plots, "mds"), "mds_free", "mds_fixed")
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

  # --- Always compute metrics_table ---
  result[["metrics_table"]] <- tryCatch(
    nm_compute_metrics(se, assay_names, condition_col),
    error = function(e) {
      warning("nm_compute_metrics() failed: ", conditionMessage(e))
      NULL
    }
  )

  # --- Always compute pc1_rank ---
  result[["pc1_rank"]] <- tryCatch(
    nm_rank_pc1(se, assay_names, condition_col, verbose = verbose),
    error = function(e) {
      warning("nm_rank_pc1() failed: ", conditionMessage(e))
      NULL
    }
  )

  # --- Always compute mds1_rank ---
  result[["mds1_rank"]] <- tryCatch(
    nm_rank_mds1(se, assay_names, condition_col, verbose = verbose),
    error = function(e) {
      warning("nm_rank_mds1() failed: ", conditionMessage(e))
      NULL
    }
  )

  non_plot <- c("metrics_table", "pc1_rank", "mds1_rank")
  n_ok   <- sum(!sapply(result[setdiff(names(result), non_plot)], is.null))
  n_fail <- length(selected) - n_ok
  if (verbose) {
    message("normalization_metrics: ", n_ok, " plot(s) generated",
            if (n_fail > 0) paste0(", ", n_fail, " failed") else ".")
  }

  # --- Export results if output_dir is set ---
  .nm_export_results(result, output_dir, export_plots, export_tables,
                     plot_width, plot_height, plot_dpi, verbose)

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


  # ---- 2. Generate all 11 quality plots at once ------------------------------

  plots <- normalization_metrics(se_nm)

  # Names of available plots + metrics_table + pc1_rank
  names(plots)  # boxplot density pcv pmad pev pca correlation mds scatter qq metrics pc1_ranking metrics_table pc1_rank


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
  plots$metrics     # group-separation metrics bar chart

  # Metrics table (data.frame, always present)
  plots$metrics_table

  # PC1 ranking (data.frame, always present)
  plots$pc1_rank                # Method, PC1_VarPct, Rank — sorted desc
  plots$pc1_ranking             # horizontal bar chart of PC1 variance


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


  # ---- 6. Standalone group-separation metrics ---------------------------------

  # Compute metrics table directly (without generating plots)
  metrics_df <- nm_compute_metrics(se_nm)
  print(metrics_df)

  # Plot metrics as a faceted bar chart
  nm_plot_metrics(se_nm)

  # PC1 variance ranking
  nm_rank_pc1(se_nm)              # data.frame with Method, PC1_VarPct, Rank
  nm_plot_pc1_ranking(se_nm)      # horizontal bar chart


  # ---- 7. Auto-export plots and tables ----------------------------------------

  # Export all plots (PNG) and tables (TSV) to a directory
  plots <- normalization_metrics(se_nm,
                                 output_dir = "./results/normalization_metrics")

  # Only tables (no plots)
  plots <- normalization_metrics(se_nm,
                                 output_dir    = "./results/normalization_metrics",
                                 export_plots  = FALSE)

  # Only plots (no tables), custom dimensions
  plots <- normalization_metrics(se_nm,
                                 output_dir    = "./results/normalization_metrics",
                                 export_tables = FALSE,
                                 plot_width    = 16,
                                 plot_height   = 10,
                                 plot_dpi      = 300)


  # ---- 8. Auto-benchmark from preprocessing (no process_proteomics needed) ----

  # Prepare SE with raw + log2 assays directly from preprocessing
  se <- nm_prepare_se(preprocessing, min_reps = 3)
  SummarizedExperiment::assayNames(se)  # "raw", "log2"

  # Auto-run all 13 normalization methods + compute metrics + plots
  plots_auto <- normalization_metrics(se, methods = "all", base_assay = "log2")

  # Or only specific methods
  plots_sub <- normalization_metrics(se,
                                      methods = c("cycloess", "vsn", "MAD",
                                                  "quantile.robust"))

  # With custom cycloess parameters
  plots_custom <- normalization_metrics(se,
                                         methods     = "all",
                                         method_args = list(
                                           cycloess = list(method = "fast",
                                                           span   = 0.8)))

  # Standalone: get only the multi-assay SE (no plots)
  se_bench <- nm_run_normalizations(se, methods = "all")
  SummarizedExperiment::assayNames(se_bench)

}
