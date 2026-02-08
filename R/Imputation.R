# =============================================================================
# Imputation Module
# =============================================================================
#
# Functions for proteomics data imputation:
#   - MNAR mask by condition
#   - Minimum value imputation (MNAR)
#   - Mixed imputation (MAR + MNAR)
#   - Pre-filtering by MNAR rules
#   - Rowname renaming from rowData
#
# Dependencies:
#   - SummarizedExperiment, S4Vectors
#
# Optional dependencies:
#   - rrcovNA (for MAR imputation with impSeqRob)
#
# Note: MNAR imputation ("min" method) does NOT require external dependencies.
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

#' MNAR imputation with minimum value
#'
#' Replaces all NA with the minimum value of the matrix.
#' Equivalent to MsCoreUtils::impute_min() without dependencies.
#'
#' @param x Numeric matrix
#' @return Matrix with NA replaced by global minimum
#' @keywords internal
.impute_min <- function(x) {
  val <- min(x, na.rm = TRUE)
  x[is.na(x)] <- val
  x
}

#' Mixed imputation: MAR/MCAR with impSeqRob, MNAR with "min"
#'
#' Two-stage imputation:
#' 1. MAR/MCAR: impSeqRob (rrcovNA) - robust and recommended
#' 2. MNAR: global minimum value (no external dependencies)
#'
#' @param x Log2 intensity matrix
#' @param condition Condition vector aligned with columns
#' @param prop_na_in_condition Proportion of NA to classify as MNAR (default: 1.0)
#' @param prop_present_in_other_condition Proportion present in other conditions (default: 0.0)
#' @param min_present_in_other_condition Minimum present values (default: 1)
#' @param require_n_other_conditions Required conditions with presence (default: 1)
#' @param mar_method MAR method: "impSeqRob" or "none" (default: "impSeqRob")
#' @param impSeqRob_args List of arguments for impSeqRob (default: list(alpha = 0.9))
#' @return List with:
#'   - x_imputed: Imputed matrix
#'   - mnar_mask: MNAR cell mask
#'   - mar_mask: MAR cell mask
#'   - summary: Imputation statistics
#' @keywords internal
.impute_mixed <- function(
    x, condition,
    prop_na_in_condition = 1.0,
    prop_present_in_other_condition = 0.0,
    min_present_in_other_condition = 1,
    require_n_other_conditions = 1,
    mar_method = c("impSeqRob", "none"),
    impSeqRob_args = list(alpha = 0.9, norm_impute = FALSE, check_data = FALSE, verbose = TRUE)
) {
  x <- as.matrix(x)
  stopifnot(ncol(x) == length(condition))
  mar_method <- match.arg(mar_method)

  # Build masks
  mnar_mask <- .mnar_mask_by_condition(
    x, condition,
    prop_na_in_condition = prop_na_in_condition,
    prop_present_in_other_condition = prop_present_in_other_condition,
    min_present_in_other_condition = min_present_in_other_condition,
    require_n_other_conditions = require_n_other_conditions
  )
  mar_mask <- is.na(x) & !mnar_mask

  # ---- Stage 1: MAR/MCAR imputation with impSeqRob ----
  x_stage1 <- x
  if (mar_method == "impSeqRob") {
    if (!requireNamespace("rrcovNA", quietly = TRUE)) {
      stop("Para mar_method = 'impSeqRob' necesitas el paquete 'rrcovNA'.")
    }
    # Validate arguments
    allowed <- c("alpha", "norm_impute", "check_data", "verbose")
    bad <- setdiff(names(impSeqRob_args), allowed)
    if (length(bad)) {
      warning("Argumentos no soportados para impSeqRob(): ",
              paste(bad, collapse = ", "), ". Se ignoran.")
    }
    args_final <- modifyList(
      list(alpha = 0.9, norm_impute = FALSE, check_data = FALSE, verbose = TRUE),
      impSeqRob_args[names(impSeqRob_args) %in% allowed]
    )

    imp1 <- do.call(rrcovNA::impSeqRob, c(list(x = x), args_final))
    x_imp1 <- if (is.list(imp1) && !is.null(imp1$x)) imp1$x else as.matrix(imp1)
    x_stage1[mar_mask] <- x_imp1[mar_mask]
  }
  # else: "none" -> leave MAR/MCAR as NA

  # ---- Stage 2: MNAR imputation with minimum value ----
  x_min_all <- .impute_min(x_stage1)

  x_final <- x_stage1
  x_final[mnar_mask] <- x_min_all[mnar_mask]

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

#' Pre-filter proteins by MNAR rules
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
#' @return List with:
#'   - keep: Logical vector of rows to keep
#'   - summary: Filtering summary
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
#' Complete imputation pipeline: MNAR pre-filtering, mixed MAR+MNAR imputation,
#' rowname renaming, and SummarizedExperiment update.
#'
#' @param se SummarizedExperiment with normalized assay (output of normalize_proteomics)
#' @param normalized_assay_name Name of the normalized assay to use (default: "normalized")
#' @param imputed_assay_name Name for the imputed assay (default: "imputed")
#' @param prop_na_mnar NA proportion threshold for MNAR classification (default: 0.51)
#' @param prop_present_mar Present proportion for MAR (default: 0.5)
#' @param min_present_mar Minimum present values for MAR (default: 1)
#' @param require_n_conditions Number of conditions required with presence (default: 1)
#' @param mar_method MAR imputation method: "impSeqRob" or "none" (default: "impSeqRob")
#' @param impSeqRob_args List of arguments for impSeqRob
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return List with:
#'   \itemize{
#'     \item se: SummarizedExperiment with imputed assay added
#'     \item x_imputed: Raw imputed matrix (for direct export, avoids SE alignment issues)
#'     \item prefilter_summary: Pre-filtering summary
#'     \item imputation_summary: Imputation statistics
#'     \item mnar_mask: MNAR mask matrix
#'     \item mar_mask: MAR mask matrix
#'   }
#'
#' @examples
#' \dontrun{
#' imp_result <- impute_proteomics(
#'   se = norm_result$se,
#'   normalized_assay_name = "normalized",
#'   imputed_assay_name = "Cycloess"
#' )
#' }
#'
#' @export
impute_proteomics <- function(
    se,
    normalized_assay_name = "normalized",
    imputed_assay_name = "imputed",
    prop_na_mnar = 0.51,
    prop_present_mar = 0.5,
    min_present_mar = 1,
    require_n_conditions = 1,
    mar_method = c("impSeqRob", "none"),
    impSeqRob_args = list(alpha = 0.9, norm_impute = FALSE,
                          check_data = FALSE, verbose = TRUE),
    verbose = TRUE
) {
  mar_method <- match.arg(mar_method)

  # Validate SE
  stopifnot(inherits(se, "SummarizedExperiment"))
  if (!normalized_assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", normalized_assay_name, "' no encontrado en SE. ",
         "Assays disponibles: ",
         paste(SummarizedExperiment::assayNames(se), collapse = ", "))
  }

  # =========================================================================
  # 1. EXTRACT NORMALIZED DATA AND CONDITIONS
  # =========================================================================

  x_norm <- SummarizedExperiment::assay(se, normalized_assay_name)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  condition_vec <- as.factor(cd$Condition)

  # =========================================================================
  # 2. PRE-FILTER BY MNAR RULES
  # =========================================================================

  if (verbose) cat("\n=== PREFILTRADO POR REGLAS MNAR ===\n")

  pf <- .prefilter_by_rules(
    x = x_norm,
    condition = condition_vec,
    prop_na_in_condition = prop_na_mnar,
    prop_present_in_other_condition = prop_present_mar,
    min_present_in_other_condition = min_present_mar,
    require_n_other_conditions = require_n_conditions
  )
  keep_rows <- pf$keep

  if (verbose) {
    cat("- Proteinas conservadas:", pf$summary$n_keep, "\n")
    cat("- Proteinas eliminadas:", pf$summary$n_drop, "\n")
  }

  x_norm_prefilt <- x_norm[keep_rows, , drop = FALSE]

  # =========================================================================
  # 3. MIXED IMPUTATION
  # =========================================================================

  if (verbose) cat("\n=== IMPUTACION MIXTA (MAR + MNAR) ===\n")

  res_impute <- .impute_mixed(
    x = x_norm_prefilt,
    condition = condition_vec,
    prop_na_in_condition = prop_na_mnar,
    prop_present_in_other_condition = prop_present_mar,
    min_present_in_other_condition = min_present_mar,
    require_n_other_conditions = require_n_conditions,
    mar_method = mar_method,
    impSeqRob_args = impSeqRob_args
  )
  x_imputed <- res_impute$x_imputed

  if (verbose) {
    cat("- NA inicial:", round(res_impute$summary$na_rate_initial * 100, 2), "%\n")
    cat("- NA despues MAR:", round(res_impute$summary$na_rate_after_mar * 100, 2), "%\n")
    cat("- NA final:", round(res_impute$summary$na_rate_final * 100, 2), "%\n")
  }

  # =========================================================================
  # 4. RENAME ROWNAMES AND ALIGN SE
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
    prefilter_summary = pf$summary,
    imputation_summary = res_impute$summary,
    mnar_mask = res_impute$mnar_mask,
    mar_mask = res_impute$mar_mask
  )
}
