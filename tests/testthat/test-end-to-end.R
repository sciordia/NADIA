# =============================================================================
# The pipeline end to end, and the modules that consume its output
# =============================================================================
#
# The dimensions asserted here are the ones the earlier phases' regression
# checks used. They are deliberately exact: a change in any of them means the
# pipeline produced a different answer, which is precisely what should never
# happen silently.

# Not skipped on CRAN/Bioconductor: this is the test that says the pipeline
# still produces the answer it is supposed to, which is exactly what a build
# machine should be checking.

test_that("the default DIA pipeline gives the expected shape", {
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    expect_equal(nrow(nadia_dia$protein_quant), 2000L)
    expect_equal(ncol(res$se_proc), 12L)
    expect_equal(nrow(res$se_proc), 1997L)   # 3 dropped by the group filter

    expect_identical(as.character(res$comparisons), c("B-A", "D-A", "D-B"))
    expect_equal(nrow(res$DEPs_results), 1997L * 3L)

    # Every protein appears exactly once per comparison.
    expect_true(all(table(res$DEPs_results$Comparison) == 1997L))
})


test_that("the missingness of the example dataset is preserved", {
    # The dataset was rebuilt with random sampling precisely so that these two
    # numbers match the full experiment. An earlier version selected the
    # best-covered proteins and had 0 % missing, which made every imputation
    # example meaningless.
    data(nadia_dia, package = "NADIA")

    quant <- as.matrix(nadia_dia$protein_quant[
        grep("^PG.Quantity_", colnames(nadia_dia$protein_quant))])
    quant[quant == 0] <- NA

    expect_gt(mean(is.na(quant)), 0.05)
    expect_lt(mean(is.na(quant)), 0.12)

    # And the differential expression must not be degenerate either.
    res <- process_proteomics(nadia_dia, verbose = FALSE)
    prop_changed <- mean(res$DEPs_results$Change != "No Change")
    expect_gt(prop_changed, 0.20)
    expect_lt(prop_changed, 0.55)
})


test_that("the documented dimensions of nadia_dia are correct", {
    # R/data.R states these; a mismatch is a documentation bug.
    data(nadia_dia, package = "NADIA")

    expect_equal(dim(nadia_dia$metadata),      c(12L,   7L))
    expect_equal(dim(nadia_dia$protein_id),    c(2000L, 52L))
    expect_equal(dim(nadia_dia$protein_quant), c(2000L, 44L))
    # R.Condition is an ordered factor, so the levels carry the display order.
    expect_s3_class(nadia_dia$metadata$R.Condition, "factor")
    expect_identical(levels(nadia_dia$metadata$R.Condition), c("A", "B", "D"))
})


test_that("the TMT and LFQ pipelines run to completion", {
    tmt <- preprocess_tmt(
        system.file("extdata", "nadia_tmt_report.tsv.gz", package = "NADIA"),
        condition_order = c("A", "B", "C", "D"), verbose = FALSE)
    lfq <- preprocess_lfq(
        system.file("extdata", "nadia_lfq_report.tsv.gz", package = "NADIA"),
        annot_path = system.file("extdata", "nadia_lfq_annotation.tsv",
                                 package = "NADIA"),
        verbose = FALSE)

    res_tmt <- process_proteomics(tmt, verbose = FALSE)
    res_lfq <- process_proteomics(lfq, verbose = FALSE)

    # 4 conditions give 6 pairwise contrasts; 2 conditions give 1.
    expect_length(res_tmt$comparisons, 6L)
    expect_length(res_lfq$comparisons, 1L)
    expect_identical(as.character(res_lfq$comparisons), "MUT-WT")

    for (r in list(res_tmt, res_lfq)) {
        expect_s3_class(r, "proteomics_result")
        expect_gt(nrow(r$DEPs_results), 0)
        expect_false(anyNA(r$DEPs_results$logFC))
    }
})


test_that("the visualisation modules accept the pipeline output unchanged", {
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    v <- volcano_highchart_list(res$DEPs_results)
    expect_length(v, length(res$comparisons))
    expect_s3_class(v[[1]], "highchart")

    b <- boxplot_highchart_list(res$BoxPlot_Input)
    expect_gt(length(b), 0)

    p <- pca_highchart_list(res$PCA_Input, modes = c("all", "any"))
    expect_gt(length(p), 0)
})

test_that("pca_highchart_list takes comparison names as modes", {
    # `modes` accepts "all", "any" and any comparison name. Before Phase 3 an
    # unrecognised mode aborted; now it warns and is skipped, which is what lets
    # a mixed vector be passed without checking each entry first.
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    p <- pca_highchart_list(res$PCA_Input, modes = c("all", "B-A"))
    expect_length(p, 2L)
    expect_s3_class(p[[1]], "highchart")

    expect_warning(q <- pca_highchart_list(res$PCA_Input, modes = "not_a_mode"),
                   "not found")
    expect_length(q, 0L)
})


test_that("the reactable modules accept the pipeline output unchanged", {
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    expect_s3_class(results_list_reactable(res$DEPs_results), "reactable")
    expect_s3_class(protein_list_reactable(nadia_dia$protein_id,
                                           metadata = nadia_dia$metadata),
                    "reactable")
    expect_s3_class(summary_list_widget(nadia_dia$metadata), "shiny.tag.list")
})

test_that("summary_list_widget degrades gracefully for TMT metadata", {
    # Proteome Discoverer does not export the Spectronaut identification counts.
    # The widget must show the columns it has instead of aborting.
    tmt <- preprocess_tmt(
        system.file("extdata", "nadia_tmt_report.tsv.gz", package = "NADIA"),
        condition_order = c("A", "B"), verbose = FALSE)

    expect_false("R.PrecursorsIdentified" %in% names(tmt$metadata))
    expect_s3_class(summary_list_widget(tmt$metadata), "shiny.tag.list")
})


test_that("the widgets ship their JavaScript, with no CDN reference", {
    # The Excel export must work offline, which means the libraries have to be
    # inside the package and served as htmlDependency objects.
    js <- list.files(system.file("js", package = "NADIA"), recursive = TRUE)
    expect_true(any(grepl("exceljs", js)))
    expect_true(any(grepl("papaparse", js)))
    expect_true(any(grepl("LICENSE", js)))   # both are MIT; the terms travel

    deps <- NADIA:::.rl_export_deps()
    expect_length(deps, 2L)
    for (d in deps) {
        expect_s3_class(d, "html_dependency")
        expect_true(file.exists(file.path(d$src$file, d$script)))
    }

    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)
    f <- file.path(tempdir(), "nadia_widget_probe.html")
    on.exit(unlink(f), add = TRUE)
    htmltools::save_html(results_list_widget(head(res$DEPs_results, 20)), f)

    html <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_false(grepl("cdnjs", html, fixed = TRUE))
    expect_false(grepl("src=\"https://", html, fixed = TRUE))
})


test_that("benchmarking scores the pipeline against the spike-in truth", {
    skip_if_not_installed("pROC")
    data(nadia_dia, package = "NADIA")

    de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results
    sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                        package = "NADIA"))
    species_df <- data.frame(Protein.IDs = sp$PG.ProteinGroups,
                             Species     = sp$PG.OrganismId)
    ev <- data.frame(Comparison = rep(c("B-A", "D-A", "D-B"), each = 2),
                     Species = c("ECOLI", "YEAST"),
                     expected_logFC = c(1, -0.58, 2, -3.3, 1, -2.72))

    b <- benchmarking_proteomics(de, ev, species_df = species_df,
                                 output_dir = NULL, verbose = FALSE)

    # The invariant that makes specificity mean anything: a prediction can never
    # move a protein between the biological classes.
    tab <- table(b$classified_df$classification, b$classified_df$truth)
    expect_equal(unname(tab["FP", "1"]), 0L)
    expect_equal(unname(tab["TN", "1"]), 0L)
    expect_equal(unname(tab["TP", "0"]), 0L)
    expect_equal(unname(tab["FN", "0"]), 0L)

    # And the pipeline must actually perform: this is a spike-in, after all.
    expect_true(all(b$metrics_table$Specificity > 0.8))
    expect_true(all(b$metrics_table$AUC > 0.7))
})
