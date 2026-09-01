# =============================================================================
# The input a functional analysis needs: differential abundance plus clusters
# =============================================================================
#
# Differential abundance says which proteins moved; Pattern Profiler says which
# pattern they follow. An enrichment or functional analysis needs both at once,
# and joining them by hand is easy to get wrong: `DEPs_results` is long by
# comparison, `long_output` is long by cluster, and merging them without saying
# what to do about that produces a silent cartesian product.
#
# The same table is available two ways. `deps_with_clusters()` builds it in
# memory, straight from the two results, with no database involved.
# `nadia_deps_with_clusters()` reads it from a `.nadia` file, where the join is
# the SQL view `v_deps_pattern_profiler` -- so a client that is not R gets the
# same table from the same definition. The two agree row for row, which the
# tests assert rather than assume.
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================


#' Rank a protein's clusters by membership
#'
#' The R side of the `pp` CTE in `v_deps_pattern_profiler`, and deliberately
#' written to mirror it: `ClusterRank` orders a protein's clusters by decreasing
#' membership, ties broken on the cluster number, and a protein that
#' `min_membership` removed from `long_output` altogether comes back from the
#' hard assignment at rank 1.
#'
#' @param pp Result of `pattern_profiler_analysis()`.
#' @return A data frame with `Protein.IDs`, `Cluster`, `Membership` and
#'   `ClusterRank`.
#' @keywords internal
#' @noRd
.deps_clusters_table <- function(pp) {
  lo <- pp$long_output

  tab <- data.frame(
    Protein.IDs = as.character(lo$FeatureID),
    Cluster     = as.integer(lo$Cluster),
    Membership  = as.numeric(lo$Membership),
    stringsAsFactors = FALSE, row.names = NULL)

  tab <- tab[order(tab$Protein.IDs, -tab$Membership, tab$Cluster), ,
             drop = FALSE]
  tab$ClusterRank <- as.integer(
    stats::ave(seq_len(nrow(tab)), tab$Protein.IDs, FUN = seq_along))

  # A protein whose highest membership fell below the threshold used for the
  # analysis has no row above, and it is not the same thing as a protein that
  # never entered the clustering. The hard assignment tells them apart.
  m <- pp$cl$membership
  if (!is.null(m) && nrow(m) > 0) {
    miss <- setdiff(rownames(m), tab$Protein.IDs)
    if (length(miss) > 0) {
      mm   <- m[miss, , drop = FALSE]
      hard <- max.col(mm, ties.method = "first")
      tab  <- rbind(tab, data.frame(
        Protein.IDs = miss,
        Cluster     = as.integer(hard),
        Membership  = as.numeric(mm[cbind(seq_along(miss), hard)]),
        ClusterRank = 1L,
        stringsAsFactors = FALSE, row.names = NULL))
    }
  }

  rownames(tab) <- NULL
  tab
}


#' Apply the filters both routes share
#'
#' Written once and called by `deps_with_clusters()` and by
#' `nadia_deps_with_clusters()`, so that the in-memory table and the one read
#' from a file cannot drift apart in what they mean.
#'
#' The order matters. The differential-abundance side is filtered first, then
#' the cluster side, and a protein whose clusters are all filtered out keeps one
#' row with no cluster rather than disappearing: the full set of proteins that
#' were tested is the background of any enrichment, and losing part of it
#' silently inflates every result computed against it.
#'
#' @param x                The joined table.
#' @param de_cols          Column names of the differential-abundance side.
#' @param assignment       `"primary"` or `"all"`.
#' @param min_membership   Extra membership threshold, or `NULL`.
#' @param comparison       Comparisons to keep, or `NULL` for all.
#' @param significant_only Keep only rows classified `Up` or `Down`.
#' @return `x`, filtered.
#' @keywords internal
#' @noRd
.deps_clusters_filter <- function(x, de_cols, assignment, min_membership,
                                  comparison, significant_only) {

  if (!is.null(comparison)) {
    known <- unique(as.character(x$Comparison))
    bad   <- setdiff(as.character(comparison), known)
    if (length(bad) > 0) {
      stop("Comparison not present in the results: ",
           paste(bad, collapse = ", "), ". Available: ",
           paste(known, collapse = ", "), ".", call. = FALSE)
    }
    x <- x[as.character(x$Comparison) %in% as.character(comparison), ,
           drop = FALSE]
  }

  if (isTRUE(significant_only)) {
    if (!"Change" %in% names(x)) {
      stop("`significant_only = TRUE` needs a 'Change' column, which the ",
           "differential-abundance results do not have.", call. = FALSE)
    }
    x <- x[!is.na(x$Change) & as.character(x$Change) != "No Change", ,
           drop = FALSE]
  }

  drop <- rep(FALSE, nrow(x))
  if (!is.null(min_membership)) {
    drop <- !is.na(x$Membership) & x$Membership < min_membership
  }
  if (identical(assignment, "primary")) {
    drop <- drop | (!is.na(x$ClusterRank) & x$ClusterRank != 1L)
  }

  # One row per protein and comparison is only a guarantee if a group that
  # loses every cluster keeps a row with none.
  key   <- paste(x$Protein.IDs, x$Comparison, sep = "\r")
  keep  <- !drop
  empty <- !(key %in% key[keep]) & !duplicated(key)

  x$Cluster[empty]     <- NA_integer_
  x$Membership[empty]  <- NA_real_
  x$ClusterRank[empty] <- NA_integer_

  x <- x[keep | empty, c(de_cols, "Cluster", "Membership", "ClusterRank"),
         drop = FALSE]
  rownames(x) <- NULL
  x
}


#' Differential abundance with the Pattern Profiler cluster of each protein
#'
#' Joins the differential-abundance results to the soft clustering, giving the
#' table a functional-enrichment or over-representation analysis needs: the
#' statistics that say which proteins moved, and the cluster that says which
#' pattern they follow.
#'
#' The join is a left join from the results, so **every protein that was tested
#' survives it**, with or without a cluster. That matters more than it looks:
#' the complete set of tested proteins is the background an enrichment is
#' computed against, and a protein with no cluster -- one that was never
#' selected for the clustering, or that the clustering dropped -- belongs in
#' that background rather than in the bin.
#'
#' Clustering is fuzzy, so a protein can belong to several clusters at once.
#' `assignment` decides what to do about it. `"primary"` keeps only the
#' dominant pattern of each protein and returns exactly one row per protein and
#' comparison, which is what most enrichment tools expect. `"all"` keeps every
#' association above the threshold, which is the honest representation of a soft
#' clustering and the one to use when the analysis downstream can weight by
#' membership.
#'
#' @param DEPs_results Differential-abundance results, the `DEPs_results`
#'   element of [process_proteomics()]. Needs at least `Protein.IDs` and
#'   `Comparison`.
#' @param pattern_profiler Result of [pattern_profiler_analysis()], or `NULL`
#'   for no clustering at all, in which case the three cluster columns are
#'   returned as `NA` and the shape of the table is unchanged.
#' @param assignment `"primary"` (default) for one row per protein and
#'   comparison, carrying the cluster of highest membership; `"all"` for one row
#'   per protein, comparison and cluster.
#' @param min_membership Optional extra membership threshold. It can only be
#'   stricter than the one the analysis ran with: pairs below that one were
#'   discarded when `long_output` was built and cannot be recovered here. A
#'   protein left with no cluster keeps its row, with `NA`.
#' @param comparison Optional character vector of comparisons to keep. The
#'   default, `NULL`, keeps them all.
#' @param significant_only Keep only the rows classified `Up` or `Down`. The
#'   classification is the one the analysis made, with its own `alpha` and
#'   `logFC_threshold`; nothing is recomputed here.
#'
#' @return A data frame with the columns of `DEPs_results` followed by
#'   `Cluster`, `Membership` and `ClusterRank` (1 is the protein's dominant
#'   cluster).
#'
#' @seealso [nadia_deps_with_clusters()] for the same table read from a
#'   `.nadia` file, and [pattern_profiler_analysis()] for the clustering.
#'
#' @examples
#' data(nadia_dia)
#' res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#' # Without a clustering the shape is the same and the clusters are NA.
#' fi <- deps_with_clusters(res$DEPs_results, NULL)
#' head(fi)
#'
#' if (requireNamespace("Mfuzz", quietly = TRUE) &&
#'     requireNamespace("Biobase", quietly = TRUE) &&
#'     requireNamespace("e1071", quietly = TRUE)) {
#'     pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
#'                                     assay_name = "Impseqrob_min",
#'                                     verbose = FALSE)
#'
#'     primary <- deps_with_clusters(res$DEPs_results, pp)
#'     nrow(primary) == nrow(res$DEPs_results)   # the background is intact
#'
#'     soft <- deps_with_clusters(res$DEPs_results, pp, assignment = "all")
#'     table(soft$ClusterRank)
#' }
#'
#' @export
deps_with_clusters <- function(DEPs_results,
                               pattern_profiler,
                               assignment       = c("primary", "all"),
                               min_membership   = NULL,
                               comparison       = NULL,
                               significant_only = FALSE) {

  assignment <- match.arg(assignment)

  if (!is.data.frame(DEPs_results)) {
    stop("`DEPs_results` must be a data frame.", call. = FALSE)
  }
  missing_cols <- setdiff(c("Protein.IDs", "Comparison"), names(DEPs_results))
  if (length(missing_cols) > 0) {
    stop("`DEPs_results` is missing the column(s): ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(pattern_profiler) &&
      !inherits(pattern_profiler, "pattern_profiler_result")) {
    stop("`pattern_profiler` must be the object returned by ",
         "pattern_profiler_analysis(), or NULL.", call. = FALSE)
  }
  if (!is.null(min_membership) &&
      (!is.numeric(min_membership) || length(min_membership) != 1 ||
       is.na(min_membership))) {
    stop("`min_membership` must be a single number, or NULL.", call. = FALSE)
  }

  de_cols <- names(DEPs_results)
  de      <- DEPs_results
  de$.row <- seq_len(nrow(de))

  if (is.null(pattern_profiler)) {
    warning("No Pattern Profiler result: the cluster columns are all NA.",
            call. = FALSE)
    cl <- data.frame(Protein.IDs = character(0), Cluster = integer(0),
                     Membership = numeric(0), ClusterRank = integer(0),
                     stringsAsFactors = FALSE)
  } else {
    cl <- .deps_clusters_table(pattern_profiler)
  }

  x <- merge(de, cl, by = "Protein.IDs", all.x = TRUE, sort = FALSE)
  x <- x[order(x$.row, x$ClusterRank), , drop = FALSE]

  .deps_clusters_filter(x, de_cols, assignment, min_membership, comparison,
                        significant_only)
}
