# =============================================================================
# deps_with_clusters() and nadia_deps_with_clusters()
# =============================================================================
#
# The join between differential abundance and the clustering exists twice: once
# in R and once as the SQL view inside the .nadia file. That is deliberate --
# a browser reading the file with DuckDB has no R to call -- but it is only
# safe while the two agree, so the first test here compares them directly.

skip_if_no_clustering <- function() {
    skip_if_not_installed("Mfuzz")
    skip_if_not_installed("Biobase")
    skip_if_not_installed("e1071")
}

skip_if_no_duck <- function() {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")
}

# Built once for the whole file: process_proteomics() plus a clustering run is
# slow enough that repeating it per test would dominate the suite.
.fi_cache <- new.env(parent = emptyenv())

fi_fixture <- function(min_membership = 0.25) {
    key <- paste0("m", min_membership)
    if (!is.null(.fi_cache[[key]])) return(.fi_cache[[key]])

    data(nadia_dia, package = "NADIA", envir = environment())
    res <- suppressMessages(process_proteomics(nadia_dia, verbose = FALSE))
    pp  <- suppressMessages(pattern_profiler_analysis(
        res$se_proc, res$DEPs_results, assay_name = "Impseqrob_min",
        auto_select_c = FALSE, c = 4, min_membership = min_membership,
        seed = 123, verbose = FALSE))

    .fi_cache[[key]] <- list(pre = nadia_dia, res = res, pp = pp)
    .fi_cache[[key]]
}

fi_file <- function(x, envir = parent.frame()) {
    f <- tempfile(fileext = ".nadia")
    suppressMessages(write_nadia(f, x$pre, x$res, pattern_profiler = x$pp,
                                 verbose = FALSE))
    withr::defer(unlink(f), envir = envir)
    f
}


test_that("the R join and the SQL view give the same table", {
    skip_if_no_clustering()
    skip_if_no_duck()

    x <- fi_fixture()
    f <- fi_file(x)

    for (mode in c("primary", "all")) {
        expect_identical(
            deps_with_clusters(x$res$DEPs_results, x$pp, assignment = mode),
            nadia_deps_with_clusters(f, assignment = mode),
            label = paste("assignment =", mode))
    }

    # And with every filter engaged at once, which is where an ordering
    # difference between the two implementations would show up.
    expect_identical(
        deps_with_clusters(x$res$DEPs_results, x$pp, assignment = "all",
                           min_membership = 0.4, comparison = "B-A",
                           significant_only = TRUE),
        nadia_deps_with_clusters(f, assignment = "all", min_membership = 0.4,
                                 comparison = "B-A", significant_only = TRUE))
})


test_that("'primary' gives one row per protein and comparison", {
    skip_if_no_clustering()

    x <- fi_fixture()
    p <- deps_with_clusters(x$res$DEPs_results, x$pp)

    expect_identical(anyDuplicated(paste(p$Protein.IDs, p$Comparison)), 0L)
    expect_true(all(is.na(p$ClusterRank) | p$ClusterRank == 1L))
})


test_that("the set of tested proteins survives the join", {
    skip_if_no_clustering()

    x <- fi_fixture()
    p <- deps_with_clusters(x$res$DEPs_results, x$pp)

    # The background of an enrichment is every protein that was tested, so a
    # left join that loses rows is a wrong background, not a tidier table.
    expect_identical(nrow(p), nrow(x$res$DEPs_results))
    expect_identical(p[, names(x$res$DEPs_results)], x$res$DEPs_results)
    expect_identical(names(p),
                     c(names(x$res$DEPs_results),
                       "Cluster", "Membership", "ClusterRank"))

    # A threshold no protein can meet empties the clusters without emptying
    # the table.
    strict <- deps_with_clusters(x$res$DEPs_results, x$pp,
                                 min_membership = 1.1)
    expect_identical(nrow(strict), nrow(x$res$DEPs_results))
    expect_true(all(is.na(strict$Cluster)))
})


test_that("'all' keeps every membership above the threshold", {
    skip_if_no_clustering()

    x <- fi_fixture()
    a <- deps_with_clusters(x$res$DEPs_results, x$pp, assignment = "all")
    p <- deps_with_clusters(x$res$DEPs_results, x$pp)

    expect_gte(nrow(a), nrow(p))

    # One comparison at a time: the clustering is global, so each comparison
    # sees the same associations.
    one <- deps_with_clusters(x$res$DEPs_results, x$pp, assignment = "all",
                              comparison = x$res$comparisons[1])
    de1 <- x$res$DEPs_results[x$res$DEPs_results$Comparison ==
                                  x$res$comparisons[1], ]
    lo  <- x$pp$long_output
    matched <- sum(lo$FeatureID %in% de1$Protein.IDs)
    unclustered <- sum(!de1$Protein.IDs %in% lo$FeatureID)
    expect_identical(nrow(one), matched + unclustered)
})


test_that("a protein below the analysis threshold keeps its dominant cluster", {
    skip_if_no_clustering()
    skip_if_no_duck()

    # min_membership is applied when long_output is built, so at 0.9 most
    # proteins vanish from it entirely. They were still clustered, and the
    # stored hard assignment is what says so.
    x <- fi_fixture(min_membership = 0.9)
    dropped <- setdiff(rownames(x$pp$cl$membership), x$pp$long_output$FeatureID)
    expect_gt(length(dropped), 0)

    p <- deps_with_clusters(x$res$DEPs_results, x$pp)
    hit <- p[p$Protein.IDs %in% dropped, ]
    expect_gt(nrow(hit), 0)
    expect_false(anyNA(hit$Cluster))

    # And the same through the file.
    f <- fi_file(x)
    expect_identical(nadia_deps_with_clusters(f), p)

    # The hard assignment agrees with the membership matrix it came from.
    m <- x$pp$cl$membership
    i <- match(hit$Protein.IDs, rownames(m))
    expect_identical(hit$Cluster, as.integer(max.col(m[i, , drop = FALSE],
                                                     ties.method = "first")))
})


test_that("a file written without pp_assignment still builds the view", {
    skip_if_no_clustering()
    skip_if_no_duck()

    x <- fi_fixture(min_membership = 0.9)
    f <- fi_file(x)

    con <- nadia_connect(f, read_only = FALSE)
    DBI::dbExecute(con, "DROP TABLE pp_assignment")
    NADIA:::.db_create_views(con)
    DBI::dbDisconnect(con, shutdown = TRUE)

    old <- nadia_deps_with_clusters(f)
    expect_identical(nrow(old), nrow(x$res$DEPs_results))

    # Without the hard assignment those proteins fall back to no cluster,
    # which is exactly the loss the table was added to prevent.
    dropped <- setdiff(rownames(x$pp$cl$membership), x$pp$long_output$FeatureID)
    expect_true(all(is.na(old$Cluster[old$Protein.IDs %in% dropped])))
})


test_that("no clustering gives the same table with no clusters", {
    x <- fi_fixture_light <- local({
        data(nadia_dia, package = "NADIA", envir = environment())
        list(pre = nadia_dia,
             res = suppressMessages(process_proteomics(nadia_dia,
                                                       verbose = FALSE)))
    })

    expect_warning(p <- deps_with_clusters(x$res$DEPs_results, NULL),
                   "all NA")
    expect_identical(nrow(p), nrow(x$res$DEPs_results))
    expect_identical(p[, names(x$res$DEPs_results)], x$res$DEPs_results)
    expect_true(all(is.na(p$Cluster)))

    skip_if_no_duck()
    f <- tempfile(fileext = ".nadia")
    suppressMessages(write_nadia(f, x$pre, x$res, verbose = FALSE))
    withr::defer(unlink(f))

    expect_warning(q <- nadia_deps_with_clusters(f), "no Pattern Profiler")
    expect_identical(q, p)
})


test_that("the arguments are validated", {
    skip_if_no_clustering()

    x <- fi_fixture()
    expect_error(deps_with_clusters(x$res$DEPs_results, x$pp,
                                    comparison = "Z-A"),
                 "not present in the results")
    expect_error(deps_with_clusters(x$res$DEPs_results, list()),
                 "pattern_profiler_analysis")
    expect_error(deps_with_clusters(x$res$DEPs_results[, "logFC", drop = FALSE],
                                    x$pp),
                 "missing the column")
    expect_error(deps_with_clusters(x$res$DEPs_results, x$pp,
                                    min_membership = "high"),
                 "single number")
})


test_that("the new view and table are visible in the file", {
    skip_if_no_clustering()
    skip_if_no_duck()

    x  <- fi_fixture()
    f  <- fi_file(x)
    tb <- nadia_tables(f)

    expect_true("pp_assignment" %in% tb$name[tb$type == "table"])
    expect_true("v_deps_pattern_profiler" %in% tb$name[tb$type == "view"])

    # The Pattern Profiler round trip is unchanged by any of this.
    expect_identical(nadia_pattern_profiler(read_nadia(f)), x$pp$long_output)

    d <- file.path(tempdir(), "fi_parquet")
    unlink(d, recursive = TRUE)
    withr::defer(unlink(d, recursive = TRUE))
    p <- suppressMessages(nadia_export_parquet(f, d, verbose = FALSE))
    expect_true("v_deps_pattern_profiler.parquet" %in% basename(p))
})


test_that("a clustering added later reaches a preprocessing-only file", {
    skip_if_no_clustering()
    skip_if_no_duck()

    # The view block used to sit after the gate that returns early when there
    # is no processing step, so this file kept the tables and lost the view.
    x <- fi_fixture()
    f <- tempfile(fileext = ".nadia")
    withr::defer(unlink(f))
    suppressMessages(write_nadia(f, x$pre, verbose = FALSE))
    suppressMessages(nadia_add_pattern_profiler(f, x$pp, verbose = FALSE))

    expect_identical(nadia_pattern_profiler(read_nadia(f)), x$pp$long_output)
})
