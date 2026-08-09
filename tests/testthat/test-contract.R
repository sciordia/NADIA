# =============================================================================
# Contract invariants: the vocabulary the rest of the package relies on
# =============================================================================
#
# These are the names, classes and levels that downstream modules assume without
# checking. Changing any of them breaks callers silently, so they are pinned.
#
# Note on the assertions: every one first requires that the column EXISTS. In a
# previous round an assertion of the form `all(x %in% y)` passed vacuously
# because `x` was a typo and evaluated to character(0), and `all(character(0))`
# is TRUE. Checking existence first is what stops that.

test_that("the four preprocessors return the same S3 contract", {
    dia <- preprocess_spectronaut(
        system.file("extdata", "nadia_dia_report.tsv.gz", package = "NADIA"),
        condition_order = c("A", "B", "D"), verbose = FALSE)
    diann <- preprocess_diann(
        system.file("extdata", "nadia_diann_report.tsv.gz", package = "NADIA"),
        condition_order = c("A", "B", "D"), verbose = FALSE)
    tmt <- preprocess_tmt(
        system.file("extdata", "nadia_tmt_report.tsv.gz", package = "NADIA"),
        condition_order = c("A", "B"), verbose = FALSE)
    lfq <- preprocess_lfq(
        system.file("extdata", "nadia_lfq_report.tsv.gz", package = "NADIA"),
        annot_path = system.file("extdata", "nadia_lfq_annotation.tsv",
                                 package = "NADIA"),
        verbose = FALSE)

    for (obj in list(dia, diann, tmt, lfq)) {
        expect_s3_class(obj, "proteomics_data")
        expect_identical(names(obj), c("metadata", "protein_id", "protein_quant"))
        expect_s3_class(obj$metadata, "data.frame")
        expect_true(all(c("R.FileName", "R.Condition", "R.Replicate") %in%
                        names(obj$metadata)))
        # Coding is the join key between the three tables and the pipeline.
        expect_true(all(paste0("PG.Quantity_", obj$metadata$Coding) %in%
                        names(obj$protein_quant)))
    }

    # The leading class differs, which is what the print methods dispatch on.
    expect_identical(class(dia)[1], "spectronaut_data")
    expect_identical(class(diann)[1], "diann_data")
    expect_identical(class(tmt)[1], "tmt_data")
    expect_identical(class(lfq)[1], "lfq_data")
})


test_that("process_proteomics returns the documented structure", {
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    expect_s3_class(res, "proteomics_result")
    expect_true(all(c("se_proc", "DEPs_results", "BoxPlot_Input", "PCA_Input",
                      "comparisons", "parameters") %in% names(res)))

    expect_s4_class(res$se_proc, "SummarizedExperiment")
    expect_identical(SummarizedExperiment::assayNames(res$se_proc),
                     c("raw", "log2", "cycloess", "Impseqrob_min"))
    expect_true(all(c("Column", "Condition", "Replicate") %in%
                    names(SummarizedExperiment::colData(res$se_proc))))
})


test_that("DEPs_results carries the documented columns", {
    data(nadia_dia, package = "NADIA")
    de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results

    required <- c("Protein.IDs", "Comparison", "Gene.Names", "logFC",
                  "P.Value", "adj.P.Val", "Change", "Assay",
                  "MissGlobal", "MissComp", "MissCND1", "MissCND2")
    expect_true(all(required %in% names(de)))
    expect_gt(nrow(de), 0)   # so the assertions below cannot pass vacuously

    expect_type(de$logFC, "double")
    expect_true(all(de$adj.P.Val >= 0 & de$adj.P.Val <= 1))
    expect_true(all(de$P.Value <= de$adj.P.Val + 1e-12))   # BH never lowers p
})


test_that("the four missingness columns mean four different things", {
    # MissGlobal covers the whole experiment, MissComp only the two conditions
    # being compared. These used to be a single column whose name promised the
    # former and delivered the latter, so readers took a comparison-level figure
    # for an experiment-level one. nadia_dia has three conditions (A, B, D), so
    # the two genuinely differ here.
    data(nadia_dia, package = "NADIA")
    de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results

    miss_cols <- c("MissGlobal", "MissComp", "MissCND1", "MissCND2")

    # They come last, in order from the widest scope to the narrowest
    expect_identical(tail(names(de), 4L), miss_cols)

    for (col in miss_cols) {
        expect_type(de[[col]], "double")
        v <- de[[col]][!is.na(de[[col]])]
        expect_true(all(v >= 0 & v <= 100))
    }

    # MissGlobal does not depend on the comparison: one value per protein
    per_protein <- tapply(de$MissGlobal, de$Protein.IDs,
                          function(x) length(unique(x)))
    expect_true(all(per_protein == 1L))

    # ... and it is not a rename of MissComp
    expect_true(any(de$MissGlobal != de$MissComp))

    # With balanced replicates MissComp is the mean of the two condition
    # percentages. The margin absorbs the rounding to two decimals, which is
    # applied to each of the three columns independently.
    expect_lt(max(abs(de$MissComp - (de$MissCND1 + de$MissCND2) / 2)), 0.01)
})


test_that("Change is a factor with exactly three levels in a fixed order", {
    # Plots, tables and colour maps index into this vocabulary by name.
    data(nadia_dia, package = "NADIA")
    de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results

    expect_true("Change" %in% names(de))
    expect_true(is.factor(de$Change))
    expect_identical(levels(de$Change), c("Up", "Down", "No Change"))
    expect_false(anyNA(de$Change))
})


test_that("BoxPlot_Input and PCA_Input carry the columns their modules require", {
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    expect_true(all(c("Column", "Assay", "Intensity", "Condition") %in%
                    names(res$BoxPlot_Input)))
    expect_gt(nrow(res$BoxPlot_Input), 0)

    expect_true(all(c("SampleID", "FeatureID", "Intensity", "Condition",
                      "sig_any") %in% names(res$PCA_Input)))
    expect_gt(nrow(res$PCA_Input), 0)

    # One adjP_ column per comparison, named as .adjp_col() builds them.
    expect_true(all(NADIA:::.adjp_col(as.character(res$comparisons)) %in%
                    names(res$PCA_Input)))
})


test_that("the imputed assay carries the method names, and holds no NAs", {
    data(nadia_dia, package = "NADIA")

    res <- process_proteomics(nadia_dia, mar_method = "Impseqrob",
                              mnar_method = "MinDet", verbose = FALSE)
    expect_true("Impseqrob_MinDet" %in%
                SummarizedExperiment::assayNames(res$se_proc))

    # The guard that stops NAs reaching limma, where they would silently reduce
    # the residual degrees of freedom.
    expect_false(anyNA(
        SummarizedExperiment::assay(res$se_proc, "Impseqrob_MinDet")))
})


test_that("imp_method = 'none' is the one case that may leave NAs", {
    data(nadia_dia, package = "NADIA")
    # normalize_proteomics() and impute_proteomics() both return a list whose
    # $se element is the SummarizedExperiment.
    norm <- normalize_proteomics(nadia_dia, verbose = FALSE)
    expect_identical(names(norm), c("se", "filter_summary", "na_overview"))

    imp <- impute_proteomics(norm$se, "cycloess", imp_method = "none",
                             verbose = FALSE)

    expect_true("none" %in% SummarizedExperiment::assayNames(imp$se))
    expect_true(anyNA(SummarizedExperiment::assay(imp$se, "none")))
})


test_that("row and column alignment is by name, never by position", {
    data(nadia_dia, package = "NADIA")
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    se <- res$se_proc
    for (a in SummarizedExperiment::assayNames(se)) {
        expect_identical(rownames(SummarizedExperiment::assay(se, a)),
                         rownames(se))
        expect_identical(colnames(SummarizedExperiment::assay(se, a)),
                         colnames(se))
    }

    # No protein is duplicated or lost between the SE and the results table.
    expect_identical(sort(unique(res$DEPs_results$Protein.IDs),
                          method = "radix"),
                     sort(rownames(se), method = "radix"))
})
