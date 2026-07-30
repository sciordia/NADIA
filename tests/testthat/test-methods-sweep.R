# =============================================================================
# Every normalisation and every imputation method, exercised end to end
# =============================================================================
#
# This is the test that would have caught the eight missing stats/utils imports
# found by R CMD check in Phase 3. Those would have broken norm_method "MAD",
# "quantile.robust" and "Rlr", imp_method "PI", the PVCA and the dispersion
# metrics -- all of them silent until the specific method was called.
#
# Skipped on CRAN because it runs every method on a real dataset.

# Not skipped on CRAN/Bioconductor. These are the tests most worth running on
# a build machine -- they are what would have caught the missing stats/utils
# imports -- and they add about 20 s to a check that runs in four minutes
# against a ten-minute limit. Optional dependencies are handled per method.

test_that("every normalisation method produces a usable assay", {
    data(nadia_dia, package = "NADIA")

    methods <- get(".NM_BENCH_METHODS", envir = asNamespace("NADIA"))
    expect_gt(length(methods), 10)   # the list must not be silently empty

    optional <- c(vsn = "vsn", Rlr = "MASS")

    for (m in methods) {
        # `[[` errors on a name that is absent, so index by position after a
        # match(); a method with no optional dependency gives NA.
        dep <- unname(optional[match(m, names(optional))])
        if (!is.na(dep) && !requireNamespace(dep, quietly = TRUE)) next

        norm <- normalize_proteomics(nadia_dia, norm_method = m, verbose = FALSE)
        se   <- norm$se

        # "log2" writes no new assay: it IS the baseline.
        assay_name <- if (m %in% c("log2", "log2Norm")) "log2" else m
        expect_true(assay_name %in% SummarizedExperiment::assayNames(se),
                    info = paste("norm_method =", m))

        mat <- SummarizedExperiment::assay(se, assay_name)
        expect_equal(dim(mat), dim(SummarizedExperiment::assay(se, "log2")),
                     info = paste("norm_method =", m))

        # Normalisation must not introduce infinities, which log2(0) would.
        expect_false(any(is.infinite(mat)), info = paste("norm_method =", m))
        # And it must not blank out everything.
        expect_true(any(is.finite(mat)), info = paste("norm_method =", m))
    }
})


test_that("every imputation method fills the matrix and leaves no NAs", {
    data(nadia_dia, package = "NADIA")
    se <- normalize_proteomics(nadia_dia, verbose = FALSE)$se

    all_methods <- get(".IMP_METHODS_ALL", envir = asNamespace("NADIA"))
    expect_gt(length(all_methods), 15)

    # Methods needing a package that may not be installed.
    optional <- c(bpca = "pcaMethods", knn = "impute", mice = "mice",
                  missForest = "missForest", Impseq = "rrcovNA",
                  Impseqrob = "rrcovNA", QRILC = "imputeLCMD",
                  MLE = "norm", MinProb = "imputeLCMD", limpa = "limpa")

    for (m in setdiff(all_methods, "none")) {
        dep <- unname(optional[match(m, names(optional))])
        if (!is.na(dep) && !requireNamespace(dep, quietly = TRUE)) next

        args <- list(se = se, normalized_assay_name = "cycloess",
                     imp_method = m, verbose = FALSE)
        if (m == "with") args$with_value <- 15

        imp <- suppressWarnings(suppressMessages(
            do.call(impute_proteomics, args)))

        new_assay <- setdiff(SummarizedExperiment::assayNames(imp$se),
                             SummarizedExperiment::assayNames(se))
        expect_length(new_assay, 1L)

        mat <- SummarizedExperiment::assay(imp$se, new_assay)
        expect_false(anyNA(mat),            info = paste("imp_method =", m))
        expect_false(any(is.infinite(mat)), info = paste("imp_method =", m))

        # The assay must always be consistent with the object that carries it.
        expect_identical(rownames(mat), rownames(imp$se),
                         info = paste("imp_method =", m))

        # With max_na_prop at its default of NULL no protein is dropped.
        expect_identical(rownames(imp$se), rownames(se),
                         info = paste("imp_method =", m))
    }
})


test_that("max_na_prop is disabled by default, for every method alike", {
    # It used not to be: `combo` reached the dispatcher without the pre-filter
    # while softHybrid and every single method went through it, so the same
    # argument meant different things depending on imp_method. The default is
    # now NULL and no method filters unless asked.
    data(nadia_dia, package = "NADIA")
    se <- normalize_proteomics(nadia_dia, verbose = FALSE)$se

    na_frac  <- rowMeans(is.na(SummarizedExperiment::assay(se, "cycloess")))
    over_lim <- sum(na_frac > 0.8)
    expect_gt(over_lim, 0)   # otherwise this test proves nothing

    for (m in c("combo", "min", "MinDet", "zero")) {
        imp <- suppressWarnings(suppressMessages(
            impute_proteomics(se, "cycloess", imp_method = m, verbose = FALSE)))
        expect_equal(nrow(imp$se), nrow(se), info = paste("imp_method =", m))
    }
})


test_that("max_na_prop still filters when given a value", {
    data(nadia_dia, package = "NADIA")
    se <- normalize_proteomics(nadia_dia, verbose = FALSE)$se

    na_frac  <- rowMeans(is.na(SummarizedExperiment::assay(se, "cycloess")))
    over_lim <- sum(na_frac > 0.8)
    expect_gt(over_lim, 0)

    imp <- suppressWarnings(suppressMessages(
        impute_proteomics(se, "cycloess", imp_method = "min",
                          max_na_prop = 0.8, verbose = FALSE)))
    expect_equal(nrow(imp$se), nrow(se) - over_lim)

    # And the proteins it removed are exactly the ones over the limit.
    dropped <- setdiff(rownames(se), rownames(imp$se))
    expect_true(all(na_frac[dropped] > 0.8))
})


test_that("MAR methods impute near the protein's own level, MNAR far below", {
    # A method landing on the wrong side of the distribution is the single most
    # consequential imputation bug, and a "does it run" test cannot see it.
    #
    # The reference is the median of the values that WERE observed in the rows
    # that have gaps -- not the median of the whole matrix. Missingness is
    # concentrated in low-abundance proteins, so the proteins being imputed sit
    # below the global median to begin with (12.8 against 14.1 here). Comparing
    # against the global median would call every MAR method "too low".
    data(nadia_dia, package = "NADIA")
    skip_if_not_installed("impute")
    skip_if_not_installed("rrcovNA")

    se     <- normalize_proteomics(nadia_dia, verbose = FALSE)$se
    x      <- SummarizedExperiment::assay(se, "cycloess")
    was_na <- is.na(x)
    gappy  <- rowSums(was_na) > 0

    reference <- median(x[gappy, ][!was_na[gappy, ]])

    filled <- function(m) {
        imp <- suppressWarnings(suppressMessages(
            impute_proteomics(se, "cycloess", imp_method = m, verbose = FALSE)))
        mat <- SummarizedExperiment::assay(imp$se, m)
        # Index with the mask restricted to the rows the method returned, so a
        # dropped row can never recycle the mask and give a meaningless answer.
        median(mat[was_na[rownames(mat), , drop = FALSE]])
    }

    # MAR: within a log2 unit of the level of the proteins being filled.
    for (m in c("knn", "Impseqrob")) {
        expect_lt(abs(filled(m) - reference), 1,
                  label = paste("MAR method", m))
    }

    # MNAR: clearly below it, by construction.
    for (m in c("min", "MinDet")) {
        expect_lt(filled(m), reference - 2, label = paste("MNAR method", m))
    }

    # And the ordering is the one the two families promise.
    expect_lt(filled("min"), filled("knn"))
})


test_that("the whole grid of normalisation x imputation completes", {
    # A smaller grid than the full 13 x 19, chosen to cover both scale groups
    # (methods that receive raw intensities and methods that receive log2) and
    # both missingness assumptions.
    data(nadia_dia, package = "NADIA")

    norms <- c("cycloess", "quantile", "MAD", "quantile.robust", "GlobalMedian")
    imps  <- c("min", "MinDet", "PI", "zero")

    for (n in norms) {
        for (i in imps) {
            res <- suppressWarnings(suppressMessages(
                process_proteomics(nadia_dia, norm_method = n, imp_method = i,
                                   verbose = FALSE)))
            label <- paste(n, i, sep = " + ")

            expect_s3_class(res, "proteomics_result")
            expect_gt(nrow(res$DEPs_results), 0)
            expect_false(anyNA(res$DEPs_results$logFC), info = label)
            expect_false(anyNA(res$DEPs_results$adj.P.Val), info = label)
            expect_identical(levels(res$DEPs_results$Change),
                             c("Up", "Down", "No Change"), info = label)
        }
    }
})
