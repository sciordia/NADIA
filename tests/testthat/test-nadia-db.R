# =============================================================================
# The .nadia format
# =============================================================================
#
# The whole point of the format is that nothing is lost and nothing is stored
# twice. Those two claims pull against each other, so the tests below are all
# variations on one question: after dropping everything derivable, does the file
# still give back exactly what went in?
#
# `identical()` rather than `all.equal()` throughout. A round trip that returns
# a double where there was an integer, or a character where there was an ordered
# factor, has lost something, even if the numbers match.

skip_if_no_duckdb <- function() {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")
}

# A file written once and reused, since writing it is the slow part.
local_nadia <- function(envir = parent.frame(), ...) {
    data(nadia_dia, package = "NADIA", envir = environment())
    res <- suppressMessages(process_proteomics(nadia_dia, verbose = FALSE, ...))
    f <- tempfile(fileext = ".nadia")
    suppressMessages(write_nadia(f, nadia_dia, res, verbose = FALSE))
    withr::defer(unlink(f), envir = envir)
    list(file = f, pre = nadia_dia, res = res)
}


test_that("the preprocessing tables survive the round trip untouched", {
    skip_if_no_duckdb()
    x  <- local_nadia()
    db <- read_nadia(x$file)

    expect_identical(db$metadata,      x$pre$metadata)
    expect_identical(db$protein_id,    x$pre$protein_id)
    expect_identical(db$protein_quant, x$pre$protein_quant)

    # protein_quant is the interesting one: its quantities are not stored as a
    # table of their own but as the `raw` assay, and its per-sample counts live
    # in a long metrics table. Getting it back identical means both pivots are
    # right, down to the column order.
    pre <- nadia_preprocessing(db)
    expect_s3_class(pre, "spectronaut_data")
    expect_identical(unclass(pre)[1:3], unclass(x$pre)[1:3])
})


test_that("the processing results survive the round trip untouched", {
    skip_if_no_duckdb()
    x  <- local_nadia()
    db <- read_nadia(x$file)

    expect_identical(db$DEPs_results,  x$res$DEPs_results)
    expect_identical(db$BoxPlot_Input, x$res$BoxPlot_Input)
    expect_identical(db$PCA_Input,     x$res$PCA_Input)
})


test_that("the SummarizedExperiment is rebuilt assay for assay", {
    skip_if_no_duckdb()
    x    <- local_nadia()
    back <- nadia_result(read_nadia(x$file))

    a <- x$res$se_proc
    b <- back$se_proc
    expect_identical(SummarizedExperiment::assayNames(b),
                     SummarizedExperiment::assayNames(a))
    for (nm in SummarizedExperiment::assayNames(a)) {
        expect_identical(SummarizedExperiment::assay(b, nm),
                         SummarizedExperiment::assay(a, nm))
    }
    expect_identical(as.data.frame(SummarizedExperiment::rowData(b)),
                     as.data.frame(SummarizedExperiment::rowData(a)))
    expect_identical(as.data.frame(SummarizedExperiment::colData(b)),
                     as.data.frame(SummarizedExperiment::colData(a)))
})


test_that("nothing derivable is stored", {
    skip_if_no_duckdb()
    x   <- local_nadia()
    con <- nadia_connect(x$file)
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

    # log2 is exactly log2(raw), so it is not in the file.
    assays <- DBI::dbGetQuery(con, "SELECT DISTINCT assay AS a FROM intensities")$a
    expect_false("log2" %in% assays)
    expect_true("raw" %in% assays)

    # Neither is the imputed assay: only the cells that were filled.
    n_imputed <- DBI::dbGetQuery(con, "SELECT count(*) AS n FROM imputed_values")$n
    n_cells   <- DBI::dbGetQuery(con,
        "SELECT count(*) AS n FROM intensities WHERE assay = 'raw'")$n
    expect_lt(n_imputed, n_cells / 2)

    # And the number of filled cells is exactly the number that were missing.
    nrm <- SummarizedExperiment::assay(x$res$se_proc,
                                       x$res$parameters$norm_method)
    expect_identical(as.integer(n_imputed), sum(is.na(nrm)))
})


test_that("imputed_values doubles as the MAR/MNAR mask", {
    skip_if_no_duckdb()
    x  <- local_nadia()
    db <- read_nadia(x$file)

    iv <- db$imputed_values
    expect_true(all(iv$branch %in% c("MAR", "MNAR", "filled")))
    expect_true(any(iv$branch == "MNAR"))

    # The masks are the ones the imputation itself produced, not a guess made
    # afterwards. Note that mnar_mask covers observed cells too -- it marks whole
    # conditions -- so only the intersection with the missing cells is a branch.
    md  <- S4Vectors::metadata(x$res$se_proc)
    nrm <- SummarizedExperiment::assay(x$res$se_proc,
                                       x$res$parameters$norm_method)
    mnar <- md$mnar_mask & is.na(nrm)
    expect_identical(sum(iv$branch == "MNAR"), sum(mnar))
})


test_that("the plotting functions work straight off the file", {
    skip_if_no_duckdb()
    x    <- local_nadia()
    db   <- read_nadia(x$file)
    back <- nadia_result(db)

    # Not just "it runs": the figure built from the file is the same object as
    # the one built in memory.
    expect_identical(volcano_highchart_list(back$DEPs_results),
                     volcano_highchart_list(x$res$DEPs_results))

    expect_s3_class(boxplot_highchart_list(back$BoxPlot_Input)[[1]], "highchart")
    expect_s3_class(pca_highchart_list(back$PCA_Input)[[1]], "highchart")
    expect_s3_class(results_list_reactable(back$DEPs_results), "reactable")
    expect_s3_class(protein_list_reactable(db$protein_id, metadata = db$metadata),
                    "reactable")
})


test_that("a file can hold the preprocessing alone", {
    skip_if_no_duckdb()
    data(nadia_dia, package = "NADIA")
    f <- tempfile(fileext = ".nadia")
    on.exit(unlink(f), add = TRUE)

    suppressMessages(write_nadia(f, nadia_dia, verbose = FALSE))
    db <- read_nadia(f)

    expect_identical(db$protein_quant, nadia_dia$protein_quant)
    expect_null(db$DEPs_results)
    expect_error(nadia_result(db), "only the preprocessing")
})


test_that("the imputed assay is stored whole when the method revises observed values", {
    skip_if_no_duckdb()
    skip_if_not_installed("limpa")
    data(nadia_dia, package = "NADIA")

    set.seed(1)
    res <- suppressMessages(process_proteomics(
        nadia_dia, imp_method = "limpa", de_method = "limpa", verbose = FALSE))
    f <- tempfile(fileext = ".nadia")
    on.exit(unlink(f), add = TRUE)
    suppressMessages(write_nadia(f, nadia_dia, res, verbose = FALSE))

    db <- read_nadia(f)
    # limpa is model-based and does revise observed cells, so storing only the
    # difference would silently lose them. The writer checks instead of assuming.
    expect_identical(unname(db$meta[["imputation_storage"]]), "full")
    expect_identical(
        SummarizedExperiment::assay(nadia_result(db)$se_proc, "limpa"),
        SummarizedExperiment::assay(res$se_proc, "limpa"))
})


test_that("norm_method = 'log2' round trips, with no normalized assay of its own", {
    skip_if_no_duckdb()
    data(nadia_dia, package = "NADIA")

    res <- suppressMessages(process_proteomics(nadia_dia, norm_method = "log2",
                                               verbose = FALSE))
    f <- tempfile(fileext = ".nadia")
    on.exit(unlink(f), add = TRUE)
    suppressMessages(write_nadia(f, nadia_dia, res, verbose = FALSE))

    # The awkward case: imputation started from log2, which is the one assay
    # never stored. Both the views and the reader have to derive it from raw.
    db <- read_nadia(f)
    expect_identical(db$BoxPlot_Input, res$BoxPlot_Input)
    expect_identical(db$PCA_Input,     res$PCA_Input)
    expect_identical(
        SummarizedExperiment::assay(nadia_result(db)$se_proc, "Impseqrob_min"),
        SummarizedExperiment::assay(res$se_proc, "Impseqrob_min"))
})


test_that("DIA-NN, TMT and LFQ round trip, with their different columns", {
    skip_if_no_duckdb()

    for (pre in list(
        preprocess_diann(
            system.file("extdata", "nadia_diann_report.tsv.gz", package = "NADIA"),
            condition_order = c("A", "B", "D"), verbose = FALSE),
        preprocess_tmt(
            system.file("extdata", "nadia_tmt_report.tsv.gz", package = "NADIA"),
            condition_order = c("A", "B"), verbose = FALSE),
        preprocess_lfq(
            system.file("extdata", "nadia_lfq_report.tsv.gz", package = "NADIA"),
            annot_path = system.file("extdata", "nadia_lfq_annotation.tsv",
                                     package = "NADIA"),
            verbose = FALSE))) {

        res <- suppressMessages(process_proteomics(pre, verbose = FALSE))
        f <- tempfile(fileext = ".nadia")
        suppressMessages(write_nadia(f, pre, res, verbose = FALSE))
        db <- read_nadia(f)

        expect_identical(db$metadata,      pre$metadata)
        expect_identical(db$protein_id,    pre$protein_id)
        expect_identical(db$protein_quant, pre$protein_quant)
        expect_identical(db$PCA_Input,     res$PCA_Input)

        # The leading class has to survive, or the object comes back as
        # something it is not and print() dispatches to the wrong method.
        expect_identical(class(nadia_preprocessing(db))[1], class(pre)[1])
        unlink(f)
    }
})

test_that("DIA-NN records its own experiment_type", {
    skip_if_no_duckdb()

    pre <- preprocess_diann(
        system.file("extdata", "nadia_diann_report.tsv.gz", package = "NADIA"),
        condition_order = c("A", "B", "D"), verbose = FALSE)
    f <- tempfile(fileext = ".nadia")
    suppressMessages(write_nadia(f, pre, verbose = FALSE))
    on.exit(unlink(f), add = TRUE)

    db <- read_nadia(f)
    expect_identical(unname(db$meta[["experiment_type"]]), "diann")
    expect_s3_class(nadia_preprocessing(db), "diann_data")
})


test_that("a Pattern Profiler run is added to an existing file", {
    skip_if_no_duckdb()
    skip_if_not_installed("Mfuzz")
    skip_if_not_installed("Biobase")
    skip_if_not_installed("e1071")

    x  <- local_nadia()
    pp <- suppressMessages(pattern_profiler_analysis(
        x$res$se_proc, x$res$DEPs_results,
        assay_name = "Impseqrob_min", verbose = FALSE))

    suppressMessages(nadia_add_pattern_profiler(x$file, pp, verbose = FALSE))
    db <- read_nadia(x$file)

    # The z-scores are stored once per protein rather than once per cluster
    # membership, so rebuilding the table exercises that split.
    expect_identical(nadia_pattern_profiler(db), pp$long_output)
    expect_s3_class(
        cluster_profile_highchart_list(nadia_pattern_profiler(db))[[1]],
        "highchart")
})


test_that("adding a Pattern Profiler run refuses to invent a file", {
    skip_if_no_duckdb()
    expect_error(
        nadia_add_pattern_profiler(tempfile(fileext = ".nadia"), list()),
        "does not exist")
})


test_that("the provenance is recorded", {
    skip_if_no_duckdb()
    x  <- local_nadia()
    db <- read_nadia(x$file)

    expect_identical(unname(db$meta[["schema_version"]]), "1")
    expect_identical(unname(db$meta[["storage_compat"]]), "v0.10.2")
    expect_identical(unname(db$meta[["experiment_type"]]), "dia")
    expect_true("NADIA" %in% db$packages$package)

    # The call, without which method_args and covariate_df cannot be recovered.
    expect_true("process" %in% db$calls$step)
    expect_match(db$calls$call[db$calls$step == "process"], "process_proteomics")

    pars <- db$parameters
    expect_true(all(c("norm_method", "imp_method", "alpha") %in%
                    pars$name[pars$step == "process"]))
})


test_that("a file from a newer schema is refused rather than misread", {
    skip_if_no_duckdb()
    x   <- local_nadia()
    con <- NADIA:::.db_open(x$file)
    DBI::dbExecute(con,
        "UPDATE nadia_meta SET value = '99' WHERE key = 'schema_version'")
    NADIA:::.db_close(con)

    expect_error(read_nadia(x$file), "schema version 99")
})


test_that("nadia_tables lists our tables and not DuckDB's catalogue", {
    skip_if_no_duckdb()
    x <- local_nadia()
    t <- nadia_tables(x$file)

    expect_true(all(c("samples", "proteins", "intensities", "de_results") %in%
                    t$name[t$type == "table"]))
    expect_true(all(c("v_metadata", "v_pca_input") %in% t$name[t$type == "view"]))
    # The system catalogue has some fifty pg_* and information_schema views.
    expect_false(any(grepl("^pg_", t$name)))
})


test_that("the tables can be exported to Parquet", {
    skip_if_no_duckdb()
    x   <- local_nadia()
    dir <- file.path(tempdir(), "nadia_pq_test")
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)

    paths <- suppressMessages(nadia_export_parquet(x$file, dir, verbose = FALSE))
    expect_true(all(file.exists(paths)))
    expect_true(any(grepl("v_pca_input", paths)))
})


test_that("writing refuses to clobber a file unless told to", {
    skip_if_no_duckdb()
    x <- local_nadia()
    expect_error(
        suppressMessages(write_nadia(x$file, x$pre, x$res, verbose = FALSE)),
        "already exists")
    expect_silent(
        suppressMessages(write_nadia(x$file, x$pre, x$res, overwrite = TRUE,
                                     verbose = FALSE)))
})


test_that("the wrong kind of object is rejected at the door", {
    skip_if_no_duckdb()
    data(nadia_dia, package = "NADIA")
    f <- tempfile(fileext = ".nadia")
    on.exit(unlink(f), add = TRUE)

    expect_error(write_nadia(f, list(a = 1)), "preprocess_spectronaut")
    expect_error(write_nadia(f, nadia_dia, result = list(a = 1)),
                 "process_proteomics")
})
