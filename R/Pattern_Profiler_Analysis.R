# =============================================================================
# Pattern Profiler Analysis: Mfuzz clustering from a SummarizedExperiment
# =============================================================================
#
# This script processes a SummarizedExperiment object together with differential
# expression results (DEPs_results) to perform soft clustering with Mfuzz.
#
# Input:
#   - se_proc: SummarizedExperiment with intensity assays
#   - DEPs_results: DataFrame with adj.P.Val per comparison
#
# Output:
#   - Pattern_Profiler_Input.parquet: table in LONG format
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------





# =============================================================================
# DATA EXTRACTION FUNCTIONS
# =============================================================================

#' Extract data from a SummarizedExperiment
#'
#' Retrieves the intensity matrix, the feature IDs and the sample metadata from a
#' SummarizedExperiment object.
#'
#' @param se_proc SummarizedExperiment object
#' @param assay_name Name of the assay to use (NULL = first available)
#'
#' @return List with:
#'   - intensity_matrix: features x samples matrix
#'   - feature_ids: vector of protein IDs
#'   - sample_metadata: DataFrame with SampleID, Condition, etc.
extract_se_data <- function(se_proc, assay_name = NULL) {


  # Validate input

  if (!inherits(se_proc, "SummarizedExperiment")) {
    stop("'se_proc' must be a SummarizedExperiment object")
  }


  # Get the assay name

  available_assays <- SummarizedExperiment::assayNames(se_proc)
  if (length(available_assays) == 0) {
    stop("The SummarizedExperiment has no assays")
  }

  if (is.null(assay_name)) {
    assay_name <- available_assays[1]
    message(sprintf("Using assay: '%s'", assay_name))
  } else if (!(assay_name %in% available_assays)) {
    stop(sprintf("Assay '%s' not found. Available: %s",
                 assay_name, paste(available_assays, collapse = ", ")))
  }

  # Extract the intensity matrix
  intensity_matrix <- SummarizedExperiment::assay(se_proc, assay_name)

  # Get the feature IDs from rowData

  row_data <- as.data.frame(SummarizedExperiment::rowData(se_proc))

  # Look for the ID column (several possible names)
  id_candidates <- c("FeatureID", "Protein.IDs", "ID", "protein_id", "ProteinID")
  id_col <- intersect(id_candidates, names(row_data))[1]

  if (!is.na(id_col)) {
    feature_ids <- as.character(row_data[[id_col]])
  } else {
    # Fall back to the rownames if there is no ID column
    feature_ids <- rownames(intensity_matrix)
    if (is.null(feature_ids)) {
      feature_ids <- paste0("Feature_", seq_len(nrow(intensity_matrix)))
    }
  }
  rownames(intensity_matrix) <- feature_ids

  # Get the sample metadata from colData
  sample_metadata <- as.data.frame(SummarizedExperiment::colData(se_proc))
  sample_metadata$SampleID <- rownames(sample_metadata)

  # Check that the Condition column exists

  if (!("Condition" %in% names(sample_metadata))) {
    # Try to find alternatives
    cond_candidates <- c("condition", "Group", "group", "Treatment", "treatment")
    cond_col <- intersect(cond_candidates, names(sample_metadata))[1]

    if (!is.na(cond_col)) {
      sample_metadata$Condition <- sample_metadata[[cond_col]]
      message(sprintf("Using '%s' as the condition column", cond_col))
    } else {
      stop("Column 'Condition' not found in colData")
    }
  }

  list(
    intensity_matrix = intensity_matrix,
    feature_ids = feature_ids,
    sample_metadata = sample_metadata,
    assay_name = assay_name
  )
}


#' Merge significance information from DEPs_results
#'
#' Converts DEPs_results (long format) into wide format with one adjP_* column
#' per comparison, and computes sig_any.
#'
#' @param feature_ids Vector of feature IDs
#' @param DEPs_results DataFrame with columns Protein.IDs, adj.P.Val, Comparison, Assay
#' @param assay_name Name of the assay used to filter DEPs_results
#' @param alpha Significance threshold (default: 0.05)
#'
#' @return DataFrame with FeatureID, adjP_*, sig_any
merge_significance_info <- function(feature_ids, DEPs_results, assay_name, alpha = 0.05) {

  DEPs_results <- as.data.frame(DEPs_results)

  # Validate the required columns
  required_cols <- c("Protein.IDs", "adj.P.Val", "Comparison", "Assay")
  missing_cols <- setdiff(required_cols, names(DEPs_results))
  if (length(missing_cols) > 0) {
    stop("Missing columns in DEPs_results: ", paste(missing_cols, collapse = ", "))
  }

  # Filter by Assay
  available_assays <- unique(DEPs_results$Assay)
  if (!(assay_name %in% available_assays)) {
    stop(sprintf("Assay '%s' not found in DEPs_results. Available: %s",
                 assay_name, paste(available_assays, collapse = ", ")))
  }

  DEPs_results <- DEPs_results[DEPs_results$Assay == assay_name, ]
  message(sprintf("   - Filtered DEPs_results by Assay = '%s' (%d rows)",
                  assay_name, nrow(DEPs_results)))

  # Pivot to wide format: one row per protein, one adjP_* column per comparison
  sig_wide <- DEPs_results %>%
    dplyr::select(Protein.IDs, Comparison, adj.P.Val) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = Comparison,
      values_from = adj.P.Val,
      names_prefix = "adjP_"
    )

  # Build the base DataFrame with all the feature_ids

  feature_info <- data.frame(
    FeatureID = feature_ids,
    stringsAsFactors = FALSE
  )

  # Merge with the significance information
  feature_info <- merge(
    feature_info,
    sig_wide,
    by.x = "FeatureID",
    by.y = "Protein.IDs",
    all.x = TRUE
  )

  # Compute sig_any (significant in any comparison)
  adjP_cols <- grep("^adjP_", names(feature_info), value = TRUE)

  if (length(adjP_cols) > 0) {
    feature_info$sig_any <- apply(
      feature_info[, adjP_cols, drop = FALSE], 1,
      function(x) any(x <= alpha, na.rm = TRUE)
    )
  } else {
    feature_info$sig_any <- FALSE
    warning("No adjP_* columns found in DEPs_results")
  }

  # Replace NA in sig_any with FALSE
  feature_info$sig_any[is.na(feature_info$sig_any)] <- FALSE

  feature_info
}


#' Filter features by significance
#'
#' @param feature_info DataFrame with significance columns
#' @param filter_mode Selection mode:
#'   - "any": features significant in AT LEAST one comparison (uses sig_any).
#'   - "all": ALL features, with no significance filtering (it does NOT mean
#'            "significant in every comparison").
#'   - "specific": significant in the comparison given by `comparison`.
#' @param comparison Specific comparison (for mode="specific")
#' @param alpha Significance threshold
#'
#' @return Vector of the selected FeatureIDs
filter_significant_features <- function(feature_info,
                                         filter_mode = c("any", "all", "specific"),
                                         comparison = NULL,
                                         alpha = 0.05) {

  filter_mode <- match.arg(filter_mode)

  if (filter_mode == "all") {
    return(feature_info$FeatureID)
  }

  if (filter_mode == "any") {
    if (!("sig_any" %in% names(feature_info))) {
      stop("Column 'sig_any' not found. Run merge_significance_info() first.")
    }
    selected <- feature_info$FeatureID[feature_info$sig_any == TRUE]
    message(sprintf("Filter 'any': %d of %d features are significant",
                    length(selected), nrow(feature_info)))
    return(selected)
  }

  # filter_mode == "specific"
  if (is.null(comparison)) {
    stop("Argument 'comparison' is required for filter_mode = 'specific'")
  }

  adjP_col <- paste0("adjP_", comparison)
  if (!(adjP_col %in% names(feature_info))) {
    available <- grep("^adjP_", names(feature_info), value = TRUE)
    stop(sprintf("Comparison '%s' not found. Available: %s",
                 comparison, paste(gsub("^adjP_", "", available), collapse = ", ")))
  }

  selected <- feature_info$FeatureID[
    !is.na(feature_info[[adjP_col]]) & feature_info[[adjP_col]] <= alpha
  ]
  message(sprintf("Specific filter '%s' (alpha=%.3f): %d features",
                  comparison, alpha, length(selected)))

  selected
}


#' Build the clustering matrix aggregated by condition
#'
#' @param intensity_matrix Intensity matrix (features x samples)
#' @param sample_metadata DataFrame with SampleID, Condition
#' @param selected_features Vector of features to include
#' @param condition_order Condition order (NULL = alphabetical order)
#' @param aggregate Aggregation method: "median" or "mean"
#'
#' @return Matrix of aggregated intensities (features x conditions)
build_clustering_matrix <- function(intensity_matrix,
                                     sample_metadata,
                                     selected_features,
                                     condition_order = NULL,
                                     aggregate = c("median", "mean")) {

  aggregate <- match.arg(aggregate)
  agg_fun <- if (aggregate == "median") median else mean

  # Keep only the selected features
  intensity_matrix <- intensity_matrix[selected_features, , drop = FALSE]

  # Get the unique conditions
  conditions <- unique(sample_metadata$Condition)

  if (is.null(condition_order)) {
    condition_order <- sort(conditions)
    message(sprintf("Condition order: %s", paste(condition_order, collapse = " -> ")))
  } else {
    # Check that every condition exists
    missing <- setdiff(condition_order, conditions)
    if (length(missing) > 0) {
      stop("Conditions not found: ", paste(missing, collapse = ", "))
    }
  }

  # Aggregate by condition
  agg_matrix <- matrix(
    NA_real_,
    nrow = length(selected_features),
    ncol = length(condition_order),
    dimnames = list(selected_features, condition_order)
  )

  for (cond in condition_order) {
    samples_in_cond <- sample_metadata$SampleID[sample_metadata$Condition == cond]
    samples_in_cond <- intersect(samples_in_cond, colnames(intensity_matrix))

    if (length(samples_in_cond) == 0) {
      warning(sprintf("No samples for condition '%s'", cond))
      next
    }

    agg_matrix[, cond] <- apply(
      intensity_matrix[, samples_in_cond, drop = FALSE], 1,
      function(x) agg_fun(x, na.rm = TRUE)
    )
  }

  agg_matrix
}


# =============================================================================
# CLUSTERING FUNCTIONS (Mfuzz)
# =============================================================================

#' Create an ExpressionSet for Mfuzz
#'
#' @param mat Expression matrix (features x conditions)
#' @param feature_info DataFrame with feature information (optional)
#'
#' @return ExpressionSet object
create_expression_set <- function(mat, feature_info = NULL) {

  # Make sure it is a matrix
  mat <- as.matrix(mat)

  # Create the AnnotatedDataFrame for the features
  if (!is.null(feature_info) && nrow(feature_info) == nrow(mat)) {
    rownames(feature_info) <- rownames(mat)
    fData <- Biobase::AnnotatedDataFrame(data = feature_info)
  } else {
    fData <- Biobase::AnnotatedDataFrame(data = data.frame(
      FeatureID = rownames(mat),
      row.names = rownames(mat)
    ))
  }

  # Create the AnnotatedDataFrame for the conditions
  pData <- Biobase::AnnotatedDataFrame(data = data.frame(
    Condition = colnames(mat),
    row.names = colnames(mat)
  ))

  # Create the ExpressionSet
  Biobase::ExpressionSet(
    assayData = mat,
    featureData = fData,
    phenoData = pData
  )
}


#' Standardise an ExpressionSet (row-wise z-score)
#'
#' @param eset ExpressionSet object
#'
#' @return Standardised ExpressionSet
standardize_eset <- function(eset) {

  # Drop rows with too many NAs
  eset_filtered <- Mfuzz::filter.NA(eset, thres = 0.25)

  n_removed <- nrow(eset) - nrow(eset_filtered)
  if (n_removed > 0) {
    message(sprintf("Removed %d features with >25%% NAs", n_removed))
  }

  # Impute the remaining NAs
  eset_filled <- Mfuzz::fill.NA(eset_filtered, mode = "knn")

  # Standardise row-wise (z-score)
  eset_std <- Mfuzz::standardise(eset_filled)

  # Discard constant features: sd = 0 yields (x - mean)/0 = NaN after
  # standardise, and those rows would break or degenerate mfuzz. They are
  # detected as non-finite rows in the standardised matrix.
  X_std <- Biobase::exprs(eset_std)
  finite_rows <- apply(X_std, 1, function(r) all(is.finite(r)))
  n_const <- sum(!finite_rows)
  if (n_const > 0) {
    message(sprintf("Removed %d constant features (sd = 0) after standardising",
                    n_const))
    eset_std <- eset_std[finite_rows, ]
  }

  eset_std
}


# -----------------------------------------------------------------------------
# Cluster evaluation metrics
# -----------------------------------------------------------------------------

#' Xie-Beni index
#'
#' Measures compactness against separation between clusters. Lower is better.
#'
#' @param X Data matrix (features x conditions).
#' @param U Membership matrix (features x clusters).
#' @param centers Centroid matrix (clusters x conditions).
#' @param m Fuzzifier exponent.
#' @return Numeric value of the index.
#' @keywords internal
#' @noRd
.xie_beni_index <- function(X, U, centers, m = 2) {

  n <- nrow(X)
  c <- nrow(centers)

  # Compactness: weighted sum of the within-cluster distances
  compactness <- 0
  for (i in seq_len(n)) {
    for (j in seq_len(c)) {
      dist_sq <- sum((X[i, ] - centers[j, ])^2)
      compactness <- compactness + (U[i, j]^m) * dist_sq
    }
  }

  # Separation: minimum distance between centroids
  min_sep <- Inf
  for (j1 in seq_len(c - 1)) {
    for (j2 in (j1 + 1):c) {
      sep <- sum((centers[j1, ] - centers[j2, ])^2)
      if (sep < min_sep) min_sep <- sep
    }
  }

  if (min_sep == 0) min_sep <- .Machine$double.eps

  xb <- compactness / (n * min_sep)
  xb
}


#' Fuzzy Partition Coefficient (FPC)
#'
#' Measures the crispness of the fuzzy partition. Higher is better.
#'
#' @param U Membership matrix (features x clusters).
#' @return Numeric value of the coefficient.
#' @keywords internal
#' @noRd
.fpc_index <- function(U) {
  n <- nrow(U)
  sum(U^2) / n
}


#' Average Maximum Membership (AMM)
#'
#' Mean of the maximum membership of each feature. Higher is better.
#'
#' @param U Membership matrix (features x clusters).
#' @return Numeric value of the mean.
#' @keywords internal
#' @noRd
.amm_index <- function(U) {
  mean(apply(U, 1, max))
}


#' Evaluate a range of cluster numbers
#'
#' @param eset_std Standardised ExpressionSet
#' @param c_range Vector of cluster numbers to evaluate
#' @param m Fuzziness parameter
#' @param seeds Seeds for reproducibility
#' @param verbose Show progress
#'
#' @return DataFrame with the metrics per number of clusters
evaluate_cluster_range <- function(eset_std,
                                    c_range = 2:10,
                                    m = NULL,
                                    seeds = c(42, 123, 456),
                                    verbose = TRUE) {

  # Estimate m if it is not supplied
  if (is.null(m)) {
    m <- Mfuzz::mestimate(eset_std)
    if (verbose) message(sprintf("Estimated m parameter: %.3f", m))
  }

  X <- Biobase::exprs(eset_std)
  results <- list()

  # One seed is set per replicate inside the loop; the user's RNG is restored on
  # exit from the function.
  old_rng <- .rng_state()
  on.exit(.rng_restore(old_rng), add = TRUE)

  for (c in c_range) {
    if (verbose) message(sprintf("Evaluating c = %d...", c))

    # NA (not 0): a seed that fails must not count as 0 in the mean, because XB
    # is minimised and a spurious 0 would bias the choice of c.
    xb_vals <- rep(NA_real_, length(seeds))
    fpc_vals <- rep(NA_real_, length(seeds))
    amm_vals <- rep(NA_real_, length(seeds))
    dmin_vals <- rep(NA_real_, length(seeds))

    for (s in seq_along(seeds)) {
      set.seed(seeds[s])

      cl <- tryCatch({
        Mfuzz::mfuzz(eset_std, c = c, m = m)
      }, error = function(e) NULL)

      if (is.null(cl)) next

      U <- cl$membership
      centers <- cl$centers

      xb_vals[s] <- .xie_beni_index(X, U, centers, m)
      fpc_vals[s] <- .fpc_index(U)
      amm_vals[s] <- .amm_index(U)

      # Dmin: minimum distance between centroids
      dists <- as.matrix(dist(centers))
      diag(dists) <- Inf
      dmin_vals[s] <- min(dists)
    }

    results[[as.character(c)]] <- data.frame(
      c = c,
      XB = mean(xb_vals, na.rm = TRUE),
      XB_sd = sd(xb_vals, na.rm = TRUE),
      FPC = mean(fpc_vals, na.rm = TRUE),
      AMM = mean(amm_vals, na.rm = TRUE),
      Dmin = mean(dmin_vals, na.rm = TRUE)
    )
  }

  do.call(rbind, results)
}


#' Select the optimal number of clusters
#'
#' @param eset_std Standardised ExpressionSet
#' @param c_range Range of cluster numbers to evaluate
#' @param m Fuzziness parameter
#' @param method Method: "xb", "consensus", "elbow"
#' @param verbose Show progress
#'
#' @return List with optimal_c, metrics, m
select_optimal_clusters <- function(eset_std,
                                     c_range = 2:10,
                                     m = NULL,
                                     method = c("xb", "consensus", "elbow"),
                                     verbose = TRUE) {

  method <- match.arg(method)

  if (is.null(m)) {
    m <- Mfuzz::mestimate(eset_std)
  }

  metrics <- evaluate_cluster_range(eset_std, c_range, m, verbose = verbose)

  if (method == "xb") {
    # Minimum Xie-Beni
    optimal_c <- metrics$c[which.min(metrics$XB)]

  } else if (method == "elbow") {
    # Elbow method on Dmin
    dmin_diff <- -diff(metrics$Dmin)
    elbow_idx <- which.max(dmin_diff) + 1
    optimal_c <- metrics$c[min(elbow_idx, nrow(metrics))]

  } else {
    # Consensus: mean of the rankings
    metrics$rank_XB <- rank(metrics$XB)
    metrics$rank_FPC <- rank(-metrics$FPC)
    metrics$rank_AMM <- rank(-metrics$AMM)
    metrics$rank_Dmin <- rank(-metrics$Dmin)
    metrics$avg_rank <- (metrics$rank_XB + metrics$rank_FPC +
                          metrics$rank_AMM + metrics$rank_Dmin) / 4
    optimal_c <- metrics$c[which.min(metrics$avg_rank)]
  }

  if (verbose) {
    message(sprintf("\nMethod '%s': optimal c = %d", method, optimal_c))
  }

  list(
    optimal_c = optimal_c,
    metrics = metrics,
    m = m
  )
}


#' Run Mfuzz clustering
#'
#' @param eset_std Standardised ExpressionSet
#' @param c Number of clusters
#' @param m Fuzziness parameter
#' @param seed Seed for reproducibility
#'
#' @return Mfuzz clustering object
run_mfuzz_clustering <- function(eset_std, c, m, seed = 42) {

  # The user's RNG is restored on exit (mfuzz depends on the seed).
  old_rng <- .rng_state()
  on.exit(.rng_restore(old_rng), add = TRUE)

  set.seed(seed)
  cl <- Mfuzz::mfuzz(eset_std, c = c, m = m)

  cl
}


# =============================================================================
# BUILDING THE OUTPUT IN LONG FORMAT
# =============================================================================

#' Build the output table in LONG format
#'
#' Creates one row per FeatureID-Cluster combination where
#' membership >= min_membership. This lets a protein appear in several clusters
#' (soft clustering).
#'
#' @param cl Mfuzz clustering object
#' @param eset_std Standardised ExpressionSet (with z-scores)
#' @param conditions Vector of condition names
#' @param min_membership Minimum membership threshold
#'
#' @return DataFrame in long format
build_long_output <- function(cl, eset_std, conditions, min_membership) {

  mem_matrix <- cl$membership  # n_proteins x n_clusters
  zscores <- Biobase::exprs(eset_std)  # n_proteins x n_conditions

  n_proteins <- nrow(mem_matrix)
  n_clusters <- ncol(mem_matrix)

  # Preallocate the list for efficiency

  long_rows <- vector("list", n_proteins * n_clusters)
  row_idx <- 0

  for (i in seq_len(n_proteins)) {
    feature_id <- rownames(mem_matrix)[i]
    memberships <- mem_matrix[i, ]

    # Find the clusters where membership >= threshold
    qualifying_clusters <- which(memberships >= min_membership)

    for (k in qualifying_clusters) {
      row_idx <- row_idx + 1

      # Build the base row
      row_data <- list(
        FeatureID = feature_id,
        Cluster = as.integer(k),
        Membership = round(memberships[k], 6)
      )

      # Add the z-scores per condition
      for (cond in conditions) {
        row_data[[cond]] <- round(zscores[i, cond], 6)
      }

      long_rows[[row_idx]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
    }
  }

  # Combine all the rows
  result <- do.call(rbind, long_rows[seq_len(row_idx)])

  # Sort by Cluster and by decreasing Membership
  result <- result[order(result$Cluster, -result$Membership), ]
  rownames(result) <- NULL

  result
}


# =============================================================================
# MAIN FUNCTION: PATTERN PROFILER ANALYSIS
# =============================================================================

#' Pattern Profiler Analysis
#'
#' Complete clustering pipeline starting from a SummarizedExperiment.
#' Writes a parquet file in LONG format for visualisation.
#'
#' @param se_proc SummarizedExperiment with intensity data
#' @param DEPs_results DataFrame with differential expression results (must have an 'Assay' column)
#' @param assay_name Name of the assay to use (default: "LoessCyc"). It is also used to filter DEPs_results.
#' @param filter_mode Filtering mode: "any" (significant in at least one
#'   comparison), "all" (ALL features, with no significance filtering),
#'   "specific" (significant in the comparison given by `comparison`)
#' @param alpha Significance threshold (default: 0.05)
#' @param comparison Specific comparison (for filter_mode="specific")
#' @param condition_order Condition order for the profiles
#' @param aggregate Aggregation method: "median" or "mean"
#' @param c_range Range of cluster numbers to evaluate
#' @param auto_select_c Automatic selection of the number of clusters (default: TRUE)
#' @param c Fixed number of clusters (if auto_select_c=FALSE)
#' @param selection_method Selection method: "xb", "consensus", "elbow"
#' @param min_membership Minimum membership threshold for inclusion in the output
#' @param output_file Path of the output parquet file. `NULL` by default, which
#'   writes nothing to disk; the long-format data is returned regardless in the
#'   `long_output` element of the result.
#' @param verbose Show progress messages
#'
#' @return List (invisible) with the clustering results: `optimal_c`, `m`,
#'   `conditions`, feature counts, `selection_metrics`, `cluster_counts`, the
#'   Mfuzz `cl` object, the standardised `eset_std` and `long_output`, the
#'   long-format data.frame that is written when `output_file` is given.
#'
#' @examples
#' \dontrun{
#' result <- pattern_profiler_analysis(
#'   se_proc = se_proc,
#'   DEPs_results = DEPs_results,
#'   filter_mode = "any",
#'   condition_order = c("A", "B", "C", "D"),
#'   min_membership = 0.25,
#'   output_file = "data-raw/Pattern_Profiler_Input.parquet"
#' )
#' }
#' @export
pattern_profiler_analysis <- function(se_proc,
                                       DEPs_results,
                                       assay_name = "LoessCyc",
                                       filter_mode = c("any", "all", "specific"),
                                       alpha = 0.05,
                                       comparison = NULL,
                                       condition_order = NULL,
                                       aggregate = c("median", "mean"),
                                       c_range = 2:10,
                                       auto_select_c = TRUE,
                                       c = NULL,
                                       selection_method = c("xb", "consensus", "elbow"),
                                       min_membership = 0.25,
                                       output_file = NULL,
                                       verbose = TRUE) {

  filter_mode <- match.arg(filter_mode)
  aggregate <- match.arg(aggregate)
  selection_method <- match.arg(selection_method)

  # Mfuzz calls exprs() and cmeans() unqualified, so it needs Biobase and e1071
  # attached on the search path. They are attached here and released on exit, so
  # that the user's session is left as it was.
  .pp_attached <- .mfuzz_deps_attach()
  on.exit(.mfuzz_deps_detach(.pp_attached), add = TRUE)

  if (verbose) message("=== Pattern Profiler Analysis ===\n")

  # -------------------------------------------------------------------------
  # 1) Extract the data from the SummarizedExperiment
  # -------------------------------------------------------------------------
  if (verbose) message("1. Extracting data from the SummarizedExperiment...")

  se_data <- extract_se_data(se_proc, assay_name)

  if (verbose) {
    message(sprintf("   - Features: %d", length(se_data$feature_ids)))
    message(sprintf("   - Samples: %d", nrow(se_data$sample_metadata)))
    message(sprintf("   - Conditions: %s",
                    paste(unique(se_data$sample_metadata$Condition), collapse = ", ")))
  }

  # -------------------------------------------------------------------------
  # 2) Merge the significance information
  # -------------------------------------------------------------------------
  if (verbose) message("\n2. Merging significance information...")

  feature_info <- merge_significance_info(
    se_data$feature_ids,
    DEPs_results,
    assay_name = se_data$assay_name,
    alpha = alpha
  )

  n_sig <- sum(feature_info$sig_any, na.rm = TRUE)
  if (verbose) {
    message(sprintf("   - Significant features (alpha=%.3f): %d", alpha, n_sig))
  }

  # -------------------------------------------------------------------------
  # 3) Filter the features
  # -------------------------------------------------------------------------
  if (verbose) message("\n3. Filtering features...")

  selected_features <- filter_significant_features(
    feature_info,
    filter_mode = filter_mode,
    comparison = comparison,
    alpha = alpha
  )

  if (length(selected_features) < 10) {
    stop("Too few features selected (< 10). Adjust the filters.")
  }

  # -------------------------------------------------------------------------
  # 4) Build the clustering matrix
  # -------------------------------------------------------------------------
  if (verbose) message("\n4. Building the clustering matrix...")

  clustering_matrix <- build_clustering_matrix(
    se_data$intensity_matrix,
    se_data$sample_metadata,
    selected_features,
    condition_order,
    aggregate
  )

  # Keep the condition order
  conditions <- colnames(clustering_matrix)

  if (verbose) {
    message(sprintf("   - Matrix: %d features x %d conditions",
                    nrow(clustering_matrix), ncol(clustering_matrix)))
  }

  # -------------------------------------------------------------------------
  # 5) Create and standardise the ExpressionSet
  # -------------------------------------------------------------------------
  if (verbose) message("\n5. Standardising the data (z-score)...")

  eset <- create_expression_set(clustering_matrix)
  eset_std <- standardize_eset(eset)

  n_final <- nrow(eset_std)
  if (verbose) {
    message(sprintf("   - Final features (after NA filtering): %d", n_final))
  }

  # -------------------------------------------------------------------------
  # 6) Select the number of clusters
  # -------------------------------------------------------------------------
  if (auto_select_c) {
    if (verbose) message("\n6. Selecting the optimal number of clusters...")

    selection <- select_optimal_clusters(
      eset_std,
      c_range = c_range,
      method = selection_method,
      verbose = verbose
    )

    optimal_c <- selection$optimal_c
    m <- selection$m
    selection_metrics <- selection$metrics

  } else {
    if (is.null(c)) {
      stop("'c' must be specified when auto_select_c = FALSE")
    }
    optimal_c <- c
    m <- Mfuzz::mestimate(eset_std)
    selection_metrics <- NULL

    if (verbose) {
      message(sprintf("\n6. Using c = %d (fixed), m = %.3f", optimal_c, m))
    }
  }

  # -------------------------------------------------------------------------
  # 7) Run the clustering
  # -------------------------------------------------------------------------
  if (verbose) message(sprintf("\n7. Running Mfuzz clustering (c=%d)...", optimal_c))

  cl <- run_mfuzz_clustering(eset_std, c = optimal_c, m = m)

  # Count the proteins per cluster (hard assignment)
  hard_assignment <- apply(cl$membership, 1, which.max)
  cluster_counts <- table(hard_assignment)

  if (verbose) {
    message("   Cluster distribution (hard assignment):")
    for (k in seq_len(optimal_c)) {
      count <- ifelse(as.character(k) %in% names(cluster_counts),
                      cluster_counts[[as.character(k)]], 0)
      message(sprintf("     Cluster %d: %d proteins", k, count))
    }
  }

  # -------------------------------------------------------------------------
  # 8) Build the output in LONG format
  # -------------------------------------------------------------------------
  if (verbose) message(sprintf("\n8. Building the LONG table (membership >= %.2f)...",
                               min_membership))

  long_output <- build_long_output(cl, eset_std, conditions, min_membership)

  n_unique_features <- length(unique(long_output$FeatureID))
  n_rows <- nrow(long_output)
  multi_cluster <- n_rows - n_unique_features

  if (verbose) {
    message(sprintf("   - Total rows: %d", n_rows))
    message(sprintf("   - Unique features: %d", n_unique_features))
    message(sprintf("   - Features in several clusters: %d", multi_cluster))
  }

  # -------------------------------------------------------------------------
  # 9) Write the parquet file (only if a path has been given)
  # -------------------------------------------------------------------------
  # With output_file = NULL nothing is written: the function must not create
  # files in the user's workspace unless the user says where. The long-format
  # data.frame is returned in the result regardless.
  if (is.null(output_file)) {
    if (verbose) message("\n9. No output_file: no file is written")
  } else {
    if (verbose) message(sprintf("\n9. Writing file: %s", output_file))

    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    arrow::write_parquet(long_output, output_file)
  }

  if (verbose) message("\n=== Analysis complete ===")

  # -------------------------------------------------------------------------
  # Return the results (invisibly)
  # -------------------------------------------------------------------------
  result <- list(
    optimal_c = optimal_c,
    m = m,
    conditions = conditions,
    n_features_input = length(selected_features),
    n_features_final = n_final,
    n_rows_output = n_rows,
    min_membership = min_membership,
    selection_metrics = selection_metrics,
    cluster_counts = as.data.frame(cluster_counts),
    output_file = output_file,
    long_output = long_output,
    cl = cl,
    eset_std = eset_std
  )

  class(result) <- c("pattern_profiler_result", "list")

  invisible(result)
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# --- Typical usage ---
# The objects 'se_proc' (SummarizedExperiment) and 'DEPs_results' (dataframe)
# are already loaded in the environment from previous pipeline steps.
#
# source("R/Pattern_Profiler_Analysis.R")
#
# # Run the analysis (uses assay 'LoessCyc' by default)
# result <- pattern_profiler_analysis(
#   se_proc = se_proc,
#   DEPs_results = DEPs_results,
#   assay_name = "LoessCyc",  # default; also filters DEPs_results by this column
#   filter_mode = "any",
#   condition_order = c("A", "B", "C", "D"),
#   c_range = 2:8,
#   selection_method = "xb",
#   min_membership = 0.25,
#   output_file = "data-raw/Pattern_Profiler_Input.parquet"
# )
#
# # Inspect the results
# result$optimal_c
# result$selection_metrics
#
# --- Using a different assay ---
# result <- pattern_profiler_analysis(
#   se_proc = se_proc,
#   DEPs_results = DEPs_results,
#   assay_name = "log2",  # use log2 instead of LoessCyc
#   filter_mode = "any",
#   condition_order = c("A", "B", "C", "D")
# )
