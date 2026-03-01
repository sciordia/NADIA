# =============================================================================
# Imputation Module
# =============================================================================
#
# Functions for proteomics data imputation:
#   - 17 imputation methods (combo + 16 individual)
#   - MNAR mask by condition (for combo mode)
#   - Pre-filtering by MNAR rules (combo) or NA proportion (single methods)
#   - Mixed combo imputation (configurable MAR + MNAR methods)
#   - Rowname renaming from rowData
#
# Dependencies:
#   - SummarizedExperiment, S4Vectors
#
# Optional dependencies (method-specific):
#   - rrcovNA  (Impseq, Impseqrob)
#   - pcaMethods (bpca) [Bioconductor]
#   - impute (knn) [Bioconductor]
#   - mice (mice)
#   - missForest (missForest)

#   - imputeLCMD (QRILC, MinProb)
#   - norm (MLE)
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# --- Null coalescing operator ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# =============================================================================
# METHOD CLASSIFICATION CONSTANTS
# =============================================================================

.IMP_METHODS_ALL <- c(
  "combo", "bpca", "knn", "mice", "missForest", "Impseq",
  "Impseqrob", "QRILC", "MLE",
  "MinDet", "MinProb", "min", "zero", "nbavg", "with", "none"
)

.IMP_METHODS_MAR <- c(
  "bpca", "knn", "mice", "missForest", "Impseq", "Impseqrob",
  "MLE", "none"
)

.IMP_METHODS_MNAR <- c(
  "QRILC", "MinDet", "MinProb", "min", "zero", "with", "none"
)

# =============================================================================
# IMPUTATION HELPERS (.imp_*)
# =============================================================================
# Each receives a numeric matrix (proteins x samples, log2) and returns
# the same matrix with NAs imputed. Additional args via method_args list.

# --- No external dependencies (6) ---

#' @keywords internal
.imp_none <- function(x, args = list()) {
  x
}

#' @keywords internal
.imp_zero <- function(x, args = list()) {
  x[is.na(x)] <- 0
  x
}

#' Minimum value imputation (global minimum)
#' @keywords internal
.imp_min <- function(x, args = list()) {
  val <- min(x, na.rm = TRUE)
  x[is.na(x)] <- val
  x
}

#' MinDet: per-column quantile imputation
#' @param args list with optional `q` (default 0.01)
#' @keywords internal
.imp_MinDet <- function(x, args = list()) {
  q <- args$q %||% 0.01
  for (j in seq_len(ncol(x))) {
    na_idx <- is.na(x[, j])
    if (any(na_idx)) {
      x[na_idx, j] <- quantile(x[, j], probs = q, na.rm = TRUE)
    }
  }
  x
}

#' nbavg: neighbor averaging within each row
#' @keywords internal
.imp_nbavg <- function(x, args = list()) {
  for (i in seq_len(nrow(x))) {
    na_idx <- which(is.na(x[i, ]))
    if (length(na_idx) == 0) next
    non_na <- which(!is.na(x[i, ]))
    if (length(non_na) == 0) next
    for (k in na_idx) {
      dists <- abs(non_na - k)
      nearest <- non_na[order(dists)][1:min(2, length(non_na))]
      x[i, k] <- mean(x[i, nearest])
    }
  }
  x
}

#' with: replace NAs with user-specified constant
#' @param args list with `with_value`
#' @keywords internal
.imp_with <- function(x, args = list()) {
  val <- args$with_value
  if (is.null(val) || is.na(val)) {
    stop("imp_method='with' requiere un valor en 'with_value'.")
  }
  x[is.na(x)] <- val
  x
}

# --- With optional dependencies (10) ---

#' bpca: Bayesian PCA imputation (pcaMethods)
#' @param args list with optional `nPcs` (default 2)
#' @keywords internal
.imp_bpca <- function(x, args = list()) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Para imp_method='bpca' necesitas 'pcaMethods'.\n",
         "  BiocManager::install('pcaMethods')")
  }
  nPcs <- args$nPcs %||% 2
  nPcs <- min(nPcs, min(dim(x)) - 1)
  res <- pcaMethods::pca(x, method = "bpca", nPcs = nPcs)
  pcaMethods::completeObs(res)
}

#' knn: k-nearest neighbors imputation (impute)
#' @param args list with optional `k` (default 10)
#' @keywords internal
.imp_knn <- function(x, args = list()) {
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop("Para imp_method='knn' necesitas 'impute'.\n",
         "  BiocManager::install('impute')")
  }
  k <- args$k %||% 10
  k <- min(k, nrow(x) - 1)
  res <- impute::impute.knn(x, k = k)
  res$data
}

#' mice: Multiple Imputation by Chained Equations
#' @param args list with optional `m` (default 5), `maxit` (default 5), `method` (default "pmm")
#' @keywords internal
.imp_mice <- function(x, args = list()) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Para imp_method='mice' necesitas 'mice'.\n",
         "  install.packages('mice')")
  }
  m      <- args$m      %||% 5
  maxit  <- args$maxit   %||% 5
  method <- args$method  %||% "pmm"
  # mice works on data.frames with rows=samples, cols=features (transpose)
  x_t <- as.data.frame(t(x))
  imp <- mice::mice(x_t, m = m, maxit = maxit, method = method, printFlag = FALSE)
  res <- t(as.matrix(mice::complete(imp, 1)))
  dimnames(res) <- dimnames(x)
  res
}

#' missForest: Random Forest imputation
#' @param args list with optional `maxiter` (default 10), `ntree` (default 100)
#' @keywords internal
.imp_missForest <- function(x, args = list()) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("Para imp_method='missForest' necesitas 'missForest'.\n",
         "  install.packages('missForest')")
  }
  maxiter <- args$maxiter %||% 10
  ntree   <- args$ntree   %||% 100
  # missForest expects rows=observations (samples), cols=variables (proteins)
  x_t <- t(x)
  res <- missForest::missForest(x_t, maxiter = maxiter, ntree = ntree, verbose = FALSE)
  out <- t(res$ximp)
  dimnames(out) <- dimnames(x)
  out
}

#' Impseq: Sequential imputation (rrcovNA)
#' @keywords internal
.imp_Impseq <- function(x, args = list()) {
  if (!requireNamespace("rrcovNA", quietly = TRUE)) {
    stop("Para imp_method='Impseq' necesitas 'rrcovNA'.\n",
         "  install.packages('rrcovNA')")
  }
  res <- rrcovNA::impSeq(x)
  as.matrix(res)
}

#' Impseqrob: Robust sequential imputation (rrcovNA)
#' @param args list with optional `alpha` (default 0.9)
#' @keywords internal
.imp_Impseqrob <- function(x, args = list()) {
  if (!requireNamespace("rrcovNA", quietly = TRUE)) {
    stop("Para imp_method='Impseqrob' necesitas 'rrcovNA'.\n",
         "  install.packages('rrcovNA')")
  }
  alpha <- args$alpha %||% 0.9
  res <- rrcovNA::impSeqRob(x, alpha = alpha, norm_impute = FALSE,
                             check_data = FALSE, verbose = FALSE)
  if (is.list(res) && !is.null(res$x)) res$x else as.matrix(res)
}

#' QRILC: Quantile Regression Imputation of Left-Censored data (imputeLCMD)
#' @param args list with optional `tune.sigma` (default 1)
#' @keywords internal
.imp_QRILC <- function(x, args = list()) {
  if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
    stop("Para imp_method='QRILC' necesitas 'imputeLCMD'.\n",
         "  install.packages('imputeLCMD')")
  }
  tune.sigma <- args$tune.sigma %||% 1
  res <- imputeLCMD::impute.QRILC(x, tune.sigma = tune.sigma)
  res[[1]]
}

#' MLE: Maximum Likelihood Estimation (norm)
#' @keywords internal
.imp_MLE <- function(x, args = list()) {
  if (!requireNamespace("norm", quietly = TRUE)) {
    stop("Para imp_method='MLE' necesitas 'norm'.\n",
         "  install.packages('norm')")
  }
  # norm works on observations(rows) x variables(cols) → transpose
  x_t <- t(x)
  s <- norm::prelim.norm(x_t)
  thetahat <- norm::em.norm(s, showits = FALSE)
  seed <- args$seed %||% 1
  norm::rngseed(seed)
  res <- norm::imp.norm(s, thetahat, x_t)
  out <- t(res)
  dimnames(out) <- dimnames(x)
  out
}

#' MinProb: Minimum Probability imputation (imputeLCMD)
#' @param args list with optional `q` (default 0.01), `tune.sigma` (default 1)
#' @keywords internal
.imp_MinProb <- function(x, args = list()) {
  if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
    stop("Para imp_method='MinProb' necesitas 'imputeLCMD'.\n",
         "  install.packages('imputeLCMD')")
  }
  q          <- args$q          %||% 0.01
  tune.sigma <- args$tune.sigma %||% 1
  imputeLCMD::impute.MinProb(x, q = q, tune.sigma = tune.sigma)
}

# =============================================================================
# DISPATCHER
# =============================================================================

#' Central dispatcher for imputation methods
#'
#' @param x Numeric matrix (proteins x samples)
#' @param method Method name (one of .IMP_METHODS_ALL minus "combo")
#' @param method_args Named list of per-method argument lists
#' @param with_value Constant value for method "with"
#' @return Imputed matrix
#' @keywords internal
.dispatch_imputation <- function(x, method, method_args = list(), with_value = NA_real_) {
  args <- method_args[[method]] %||% list()
  if (method == "with") args$with_value <- with_value

  switch(method,
    "none"       = .imp_none(x, args),
    "zero"       = .imp_zero(x, args),
    "min"        = .imp_min(x, args),
    "MinDet"     = .imp_MinDet(x, args),
    "nbavg"      = .imp_nbavg(x, args),
    "with"       = .imp_with(x, args),
    "bpca"       = .imp_bpca(x, args),
    "knn"        = .imp_knn(x, args),
    "mice"       = .imp_mice(x, args),
    "missForest" = .imp_missForest(x, args),
    "Impseq"     = .imp_Impseq(x, args),
    "Impseqrob"  = .imp_Impseqrob(x, args),
"QRILC"      = .imp_QRILC(x, args),
    "MLE"        = .imp_MLE(x, args),
    "MinProb"    = .imp_MinProb(x, args),
    stop("Metodo de imputacion desconocido: '", method, "'. ",
         "Metodos disponibles: ", paste(.IMP_METHODS_ALL, collapse = ", "))
  )
}

# =============================================================================
# INTERNAL FUNCTIONS
# =============================================================================

#' MNAR mask by condition with evidence in other conditions
#'
#' Identifies MNAR candidate cells: NA in one condition but with
#' sufficient presence in other conditions.
#'
#' @param x Intensity matrix (proteins x samples)
#' @param condition Condition vector aligned with columns
#' @param prop_na_in_condition Minimum proportion of NA in the condition (default: 1.0 = 100%)
#' @param prop_present_in_other_condition Minimum proportion of non-NA in other condition (default: 0.0)
#' @param min_present_in_other_condition Minimum absolute non-NA count in other condition (default: 1)
#' @param require_n_other_conditions Number of other conditions that must meet criteria (default: 1)
#' @param drop_empty_levels Drop empty factor levels (default: TRUE)
#' @return Logical matrix of same dimensions indicating MNAR cells
#' @keywords internal
.mnar_mask_by_condition <- function(
    x, condition,
    prop_na_in_condition = 1.0,
    prop_present_in_other_condition = 0.0,
    min_present_in_other_condition = 1,
    require_n_other_conditions = 1,
    drop_empty_levels = TRUE
) {
  stopifnot(is.matrix(x) || is.data.frame(x))
  x <- as.matrix(x)
  stopifnot(ncol(x) == length(condition))
  stopifnot(prop_na_in_condition >= 0 && prop_na_in_condition <= 1)
  stopifnot(prop_present_in_other_condition >= 0 && prop_present_in_other_condition <= 1)
  stopifnot(min_present_in_other_condition >= 0)
  stopifnot(require_n_other_conditions >= 1)

  condition <- as.factor(condition)
  if (drop_empty_levels) condition <- droplevels(condition)

  mnar <- matrix(FALSE, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))

  for (g in levels(condition)) {
    jg <- which(condition == g)
    jn <- which(condition != g)
    if (length(jg) == 0) next

    # % NA within g
    frac_na_g <- rowMeans(is.na(x[, jg, drop = FALSE]))
    cond_na_ok <- frac_na_g >= prop_na_in_condition

    # Evidence in other conditions (evaluated per condition)
    present_ok_n <- integer(nrow(x))
    for (h in setdiff(levels(condition), g)) {
      jh <- which(condition == h)
      if (length(jh) == 0) next
      frac_present_h <- rowMeans(!is.na(x[, jh, drop = FALSE]))
      count_present_h <- rowSums(!is.na(x[, jh, drop = FALSE]))
      ok_h <- (frac_present_h >= prop_present_in_other_condition) &
        (count_present_h >= min_present_in_other_condition)
      present_ok_n <- present_ok_n + as.integer(ok_h)
    }
    cond_present_ok <- present_ok_n >= require_n_other_conditions

    rows_mnar_g <- cond_na_ok & cond_present_ok
    if (any(rows_mnar_g)) mnar[rows_mnar_g, jg] <- TRUE
  }
  mnar
}

#' Mixed combo imputation: configurable MAR + MNAR methods
#'
#' Two-stage imputation:
#' 1. MAR/MCAR: dispatched to mar_method
#' 2. MNAR: dispatched to mnar_method
#'
#' @param x Log2 intensity matrix
#' @param condition Condition vector aligned with columns
#' @param prop_na_in_condition Proportion of NA to classify as MNAR (default: 1.0)
#' @param prop_present_in_other_condition Proportion present in other conditions (default: 0.0)
#' @param min_present_in_other_condition Minimum present values (default: 1)
#' @param require_n_other_conditions Required conditions with presence (default: 1)
#' @param mar_method MAR imputation method (default: "Impseqrob")
#' @param mnar_method MNAR imputation method (default: "min")
#' @param method_args Named list of per-method argument lists
#' @param with_value Constant for method "with"
#' @return List with x_imputed, mnar_mask, mar_mask, summary
#' @keywords internal
.impute_combo <- function(
    x, condition,
    prop_na_in_condition = 1.0,
    prop_present_in_other_condition = 0.0,
    min_present_in_other_condition = 1,
    require_n_other_conditions = 1,
    mar_method = "Impseqrob",
    mnar_method = "min",
    method_args = list(),
    with_value = NA_real_
) {
  x <- as.matrix(x)
  stopifnot(ncol(x) == length(condition))

  # Validate methods
  if (!mar_method %in% c(.IMP_METHODS_MAR, .IMP_METHODS_MNAR)) {
    stop("mar_method '", mar_method, "' no reconocido. Opciones: ",
         paste(c(.IMP_METHODS_MAR, .IMP_METHODS_MNAR), collapse = ", "))
  }
  if (!mnar_method %in% c(.IMP_METHODS_MNAR, .IMP_METHODS_MAR)) {
    stop("mnar_method '", mnar_method, "' no reconocido. Opciones: ",
         paste(c(.IMP_METHODS_MNAR, .IMP_METHODS_MAR), collapse = ", "))
  }

  # Build masks
  mnar_mask <- .mnar_mask_by_condition(
    x, condition,
    prop_na_in_condition = prop_na_in_condition,
    prop_present_in_other_condition = prop_present_in_other_condition,
    min_present_in_other_condition = min_present_in_other_condition,
    require_n_other_conditions = require_n_other_conditions
  )
  mar_mask <- is.na(x) & !mnar_mask

  # ---- Stage 1: MAR imputation ----
  x_stage1 <- x
  if (mar_method != "none" && any(mar_mask)) {
    x_imp_mar <- .dispatch_imputation(x, mar_method, method_args, with_value)
    x_stage1[mar_mask] <- x_imp_mar[mar_mask]
  }

  # ---- Stage 2: MNAR imputation ----
  x_final <- x_stage1
  if (mnar_method != "none" && any(mnar_mask)) {
    x_imp_mnar <- .dispatch_imputation(x_stage1, mnar_method, method_args, with_value)
    x_final[mnar_mask] <- x_imp_mnar[mnar_mask]
  }

  # NA summary
  na0 <- mean(is.na(x))
  na1 <- mean(is.na(x_stage1))
  naF <- mean(is.na(x_final))

  list(
    x_imputed = x_final,
    mnar_mask = mnar_mask,
    mar_mask = mar_mask,
    summary = list(
      na_rate_initial = na0,
      na_rate_after_mar = na1,
      na_rate_final = naF,
      pct_na_marked_mnar = ifelse(any(is.na(x)), mean(mnar_mask[is.na(x)]), NA_real_)
    )
  )
}

#' Pre-filter proteins by MNAR rules (for combo mode)
#'
#' Keeps a protein if:
#' (A) It has MNAR in >= 1 condition, OR
#' (B) It has sufficient signal in >= require_n_other_conditions conditions
#'
#' @param x Log2 intensity matrix
#' @param condition Condition vector
#' @param prop_na_in_condition NA proportion for MNAR (default: 0.51)
#' @param prop_present_in_other_condition Required present proportion (default: 0.5)
#' @param min_present_in_other_condition Minimum present values (default: 1)
#' @param require_n_other_conditions Required conditions (default: 1)
#' @return List with keep, summary
#' @keywords internal
.prefilter_by_rules <- function(
    x, condition,
    prop_na_in_condition = 0.51,
    prop_present_in_other_condition = 0.5,
    min_present_in_other_condition = 1,
    require_n_other_conditions = 1
) {
  cond <- as.factor(condition)

  # 1) MNAR mask
  mnar_mask <- .mnar_mask_by_condition(
    x, cond,
    prop_na_in_condition = prop_na_in_condition,
    prop_present_in_other_condition = prop_present_in_other_condition,
    min_present_in_other_condition = min_present_in_other_condition,
    require_n_other_conditions = require_n_other_conditions
  )
  row_has_mnar <- rowSums(mnar_mask) > 0

  # 2) Presence by condition
  levs <- levels(cond)
  present_ok_mat <- sapply(levs, function(g) {
    jg <- which(cond == g)
    if (length(jg) == 0) return(rep(FALSE, nrow(x)))
    frac_present_g <- rowMeans(!is.na(x[, jg, drop = FALSE]))
    count_present_g <- rowSums(!is.na(x[, jg, drop = FALSE]))
    (frac_present_g >= prop_present_in_other_condition) &
      (count_present_g >= min_present_in_other_condition)
  })
  if (!is.matrix(present_ok_mat)) present_ok_mat <- cbind(present_ok_mat)

  n_conditions_with_presence <- rowSums(present_ok_mat)
  row_has_presence <- n_conditions_with_presence >= require_n_other_conditions

  # 3) Final vector
  keep <- row_has_mnar | row_has_presence

  list(
    keep = keep,
    row_has_mnar = row_has_mnar,
    n_conditions_with_presence = n_conditions_with_presence,
    present_ok_by_condition = `colnames<-`(present_ok_mat, levs),
    mnar_mask = mnar_mask,
    summary = list(
      n_total = nrow(x),
      n_keep = sum(keep),
      n_drop = sum(!keep),
      n_with_MNAR = sum(row_has_mnar),
      n_with_presence_rule = sum(row_has_presence)
    )
  )
}

#' Pre-filter proteins by NA proportion (for single methods)
#'
#' Removes proteins with more than max_na_prop fraction of NAs.
#'
#' @param x Numeric matrix
#' @param max_na_prop Maximum NA proportion (default: 0.8)
#' @return List with keep (logical vector), summary
#' @keywords internal
.prefilter_by_na_prop <- function(x, max_na_prop = 0.8) {
  na_frac <- rowMeans(is.na(x))
  keep <- na_frac <= max_na_prop

  list(
    keep = keep,
    summary = list(
      n_total = nrow(x),
      n_keep = sum(keep),
      n_drop = sum(!keep),
      max_na_prop = max_na_prop
    )
  )
}

#' Rename rownames using IDs column
#'
#' @param x_df Matrix or data frame with rownames as ProteinGroups
#' @param rd rowData data frame with IDs column
#' @param id_col Name of IDs column (default: "IDs")
#' @return Data frame with renamed rownames
#' @keywords internal
.rename_rownames_from_rd <- function(x_df, rd, id_col = "IDs") {
  stopifnot(all(rownames(x_df) %in% rownames(rd)))
  ids <- rd[rownames(x_df), id_col, drop = TRUE] |> as.character()
  if (anyNA(ids)) {
    pg_missing <- rownames(x_df)[is.na(ids)]
    stop(sprintf("IDs faltantes para %d ProteinGroups (ejemplos: %s)",
                 length(pg_missing), paste(head(pg_missing, 10), collapse = ", ")))
  }
  if (anyDuplicated(ids)) {
    warning(sprintf("IDs duplicados detectados: %d. Se aplicara make.unique().",
                    sum(duplicated(ids))))
    ids <- make.unique(ids)
  }
  x_out <- as.data.frame(x_df, check.names = FALSE)
  rownames(x_out) <- ids
  x_out
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Impute proteomics data
#'
#' Complete imputation pipeline with 17 methods. Three pathways:
#'
#' - `imp_method = "none"`: no imputation
#' - `imp_method = "combo"`: two-stage MAR+MNAR (configurable mar_method/mnar_method)
#' - Any other method: apply to all NAs (no MAR/MNAR distinction)
#'
#' @param se SummarizedExperiment with normalized assay
#' @param normalized_assay_name Name of the normalized assay to use (default: "cycloess")
#' @param imputed_assay_name Name for the imputed assay. If NULL, auto-generated:
#'   combo -> "{mar_method}_{mnar_method}", single -> "{imp_method}"
#' @param imp_method Imputation method (default: "combo"). One of:
#'   "combo", "bpca", "knn", "mice", "missForest", "Impseq", "Impseqrob",
#'   "QRILC", "MLE", "MinDet", "MinProb", "min", "zero",
#'   "nbavg", "with", "none"
#' @param mar_method MAR method for combo mode (default: "Impseqrob")
#' @param mnar_method MNAR method for combo mode (default: "min")
#' @param prop_na_mnar NA proportion threshold for MNAR classification (default: 0.51)
#' @param prop_present_mar Present proportion for MAR (default: 0.5)
#' @param min_present_mar Minimum present values for MAR (default: 1)
#' @param require_n_conditions Number of conditions required with presence (default: 1)
#' @param max_na_prop Maximum NA proportion for single-method pre-filtering (default: 0.8)
#' @param method_args Named list of per-method argument lists (e.g.,
#'   list(knn = list(k = 10), bpca = list(nPcs = 3)))
#' @param with_value Constant value for imp_method="with" (default: NA_real_)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return List with:
#'   \itemize{
#'     \item se: SummarizedExperiment with imputed assay added
#'     \item x_imputed: Raw imputed matrix (for direct export)
#'     \item prefilter_summary: Pre-filtering summary
#'     \item imputation_summary: Imputation statistics
#'     \item mnar_mask: MNAR mask matrix (combo only, NULL otherwise)
#'     \item mar_mask: MAR mask matrix (combo only, NULL otherwise)
#'   }
#'
#' @examples
#' \dontrun{
#' # Combo mode (default, backward compatible)
#' imp_result <- impute_proteomics(
#'   se = norm_result$se,
#'   normalized_assay_name = "cycloess"
#' )
#'
#' # Single method
#' imp_result <- impute_proteomics(
#'   se = norm_result$se,
#'   imp_method = "knn",
#'   method_args = list(knn = list(k = 15))
#' )
#'
#' # Combo with custom MAR/MNAR
#' imp_result <- impute_proteomics(
#'   se = norm_result$se,
#'   imp_method = "combo",
#'   mar_method = "knn",
#'   mnar_method = "MinProb"
#' )
#' }
#'
#' @export
impute_proteomics <- function(
    se,
    normalized_assay_name = "cycloess",
    imputed_assay_name    = NULL,
    imp_method            = "combo",
    mar_method            = "Impseqrob",
    mnar_method           = "min",
    prop_na_mnar          = 0.51,
    prop_present_mar      = 0.5,
    min_present_mar       = 1,
    require_n_conditions  = 1,
    max_na_prop           = 0.8,
    method_args           = list(),
    with_value            = NA_real_,
    verbose               = TRUE
) {
  # Validate SE
  stopifnot(inherits(se, "SummarizedExperiment"))
  if (!normalized_assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", normalized_assay_name, "' no encontrado en SE. ",
         "Assays disponibles: ",
         paste(SummarizedExperiment::assayNames(se), collapse = ", "))
  }

  # Validate imp_method
  if (!imp_method %in% .IMP_METHODS_ALL) {
    stop("imp_method '", imp_method, "' no reconocido. Opciones: ",
         paste(.IMP_METHODS_ALL, collapse = ", "))
  }

  # Auto-generate assay name
  if (is.null(imputed_assay_name)) {
    imputed_assay_name <- if (imp_method == "combo") {
      paste0(mar_method, "_", mnar_method)
    } else {
      imp_method
    }
  }

  # =========================================================================
  # 1. EXTRACT NORMALIZED DATA AND CONDITIONS
  # =========================================================================

  x_norm <- SummarizedExperiment::assay(se, normalized_assay_name)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  condition_vec <- as.factor(cd$Condition)

  # =========================================================================
  # PATH 1: NO IMPUTATION
  # =========================================================================

  if (imp_method == "none") {
    if (verbose) cat("\n=== IMPUTACION: none (sin imputar) ===\n")

    x_imputed <- x_norm
    pf_summary <- list(n_total = nrow(x_norm), n_keep = nrow(x_norm), n_drop = 0)
    imp_summary <- list(
      na_rate_initial = mean(is.na(x_norm)),
      na_rate_final = mean(is.na(x_norm))
    )
    mnar_mask <- NULL
    mar_mask  <- NULL

  # =========================================================================
  # PATH 2: COMBO (MAR + MNAR)
  # =========================================================================

  } else if (imp_method == "combo") {
    if (verbose) cat("\n=== PREFILTRADO POR REGLAS MNAR ===\n")

    pf <- .prefilter_by_rules(
      x = x_norm,
      condition = condition_vec,
      prop_na_in_condition = prop_na_mnar,
      prop_present_in_other_condition = prop_present_mar,
      min_present_in_other_condition = min_present_mar,
      require_n_other_conditions = require_n_conditions
    )

    if (verbose) {
      cat("- Proteinas conservadas:", pf$summary$n_keep, "\n")
      cat("- Proteinas eliminadas:", pf$summary$n_drop, "\n")
    }

    x_norm_prefilt <- x_norm[pf$keep, , drop = FALSE]

    if (verbose) cat("\n=== IMPUTACION COMBO (MAR:", mar_method, "+ MNAR:", mnar_method, ") ===\n")

    res_impute <- .impute_combo(
      x = x_norm_prefilt,
      condition = condition_vec,
      prop_na_in_condition = prop_na_mnar,
      prop_present_in_other_condition = prop_present_mar,
      min_present_in_other_condition = min_present_mar,
      require_n_other_conditions = require_n_conditions,
      mar_method = mar_method,
      mnar_method = mnar_method,
      method_args = method_args,
      with_value = with_value
    )

    x_imputed   <- res_impute$x_imputed
    pf_summary  <- pf$summary
    imp_summary <- res_impute$summary
    mnar_mask   <- res_impute$mnar_mask
    mar_mask    <- res_impute$mar_mask

    if (verbose) {
      cat("- NA inicial:", round(imp_summary$na_rate_initial * 100, 2), "%\n")
      cat("- NA despues MAR:", round(imp_summary$na_rate_after_mar * 100, 2), "%\n")
      cat("- NA final:", round(imp_summary$na_rate_final * 100, 2), "%\n")
    }

  # =========================================================================
  # PATH 3: SINGLE METHOD (all NAs)
  # =========================================================================

  } else {
    if (verbose) cat("\n=== PREFILTRADO POR PROPORCION NA (max:", max_na_prop, ") ===\n")

    pf <- .prefilter_by_na_prop(x_norm, max_na_prop = max_na_prop)

    if (verbose) {
      cat("- Proteinas conservadas:", pf$summary$n_keep, "\n")
      cat("- Proteinas eliminadas:", pf$summary$n_drop, "\n")
    }

    x_norm_prefilt <- x_norm[pf$keep, , drop = FALSE]

    if (verbose) cat("\n=== IMPUTACION:", imp_method, "===\n")

    na_before <- mean(is.na(x_norm_prefilt))
    x_imputed <- .dispatch_imputation(x_norm_prefilt, imp_method, method_args, with_value)
    na_after <- mean(is.na(x_imputed))

    pf_summary  <- pf$summary
    imp_summary <- list(na_rate_initial = na_before, na_rate_final = na_after)
    mnar_mask   <- NULL
    mar_mask    <- NULL

    if (verbose) {
      cat("- NA inicial:", round(na_before * 100, 2), "%\n")
      cat("- NA final:", round(na_after * 100, 2), "%\n")
    }
  }

  # =========================================================================
  # RENAME ROWNAMES AND ALIGN SE
  # =========================================================================

  if (verbose) cat("\n=== ACTUALIZANDO SUMMARIZEDEXPERIMENT ===\n")

  rd <- as.data.frame(SummarizedExperiment::rowData(se))

  # Rename rownames of x_imputed using IDs
  if ("IDs" %in% names(rd)) {
    rd_subset <- rd[rownames(x_imputed), , drop = FALSE]
    x_imputed_ids <- .rename_rownames_from_rd(x_imputed, rd_subset, id_col = "IDs")
  } else {
    x_imputed_ids <- x_imputed
  }

  # Align SE with imputed matrix
  common_ids <- intersect(rownames(se), rownames(x_imputed_ids))
  if (length(common_ids) == 0) {
    # Fallback: try matching via Protein.IDs column
    common_ids <- intersect(rd$Protein.IDs, rownames(x_imputed))
    if (length(common_ids) > 0) {
      # Use match() to get first occurrence only (avoids duplicate expansion with %in%)
      idx <- match(common_ids, rd$Protein.IDs)
      se_subset <- se[idx, ]
      mat <- x_imputed[common_ids, colnames(se_subset), drop = FALSE]
    } else {
      stop("No hay IDs en comun entre SE y matriz imputada")
    }
  } else {
    se_subset <- se[common_ids, ]
    mat <- as.matrix(x_imputed_ids[common_ids, colnames(se_subset), drop = FALSE])
  }

  # Assertion: SE subset must not exceed imputed matrix
  if (nrow(se_subset) > nrow(x_imputed)) {
    warning("SE alineado tiene ", nrow(se_subset), " filas vs ",
            nrow(x_imputed), " en matriz imputada. Ajustando.")
    se_subset <- se[rownames(x_imputed), ]
    mat <- as.matrix(x_imputed[, colnames(se_subset), drop = FALSE])
  }

  storage.mode(mat) <- "double"

  # Add imputed assay with the provided name
  SummarizedExperiment::assay(se_subset, imputed_assay_name) <- mat

  if (verbose) {
    cat("- Assays disponibles:",
        paste(SummarizedExperiment::assayNames(se_subset), collapse = ", "), "\n")
    cat("- Proteinas finales:", nrow(se_subset), "\n")
  }

  # =========================================================================
  # RETURN
  # =========================================================================

  list(
    se = se_subset,
    x_imputed = x_imputed,
    prefilter_summary = pf_summary,
    imputation_summary = imp_summary,
    mnar_mask = mnar_mask,
    mar_mask = mar_mask
  )
}
