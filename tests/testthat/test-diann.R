# =============================================================================
# preprocess_diann(): the DIA-NN protein-group matrix
# =============================================================================
#
# DIA-NN is the fourth input format and the one that carries the least: a
# protein-group matrix holds intensities and six annotation columns, and nothing
# per sample. These tests pin down what that means downstream -- which metric
# families are fabricated, which are left out, and that neither choice breaks
# the pipeline or the round trip through a .nadia file.

.diann_report <- function() {
    system.file("extdata", "nadia_diann_report.tsv.gz", package = "NADIA")
}

test_that("the design is read from the Abundance: column names", {
    d <- preprocess_diann(.diann_report(), condition_order = c("A", "B", "D"),
                          verbose = FALSE)

    expect_identical(d$metadata$Coding,
                     c(paste0("A_", 1:4), paste0("B_", 1:4), paste0("D_", 1:4)))
    expect_s3_class(d$metadata$R.Condition, "ordered")
    expect_identical(levels(d$metadata$R.Condition), c("A", "B", "D"))
    expect_type(d$metadata$R.Replicate, "integer")

    # Coding must reproduce the PG.Quantity_ suffixes exactly: that join is
    # what .prepare_metadata()/.prepare_protein_data() rely on.
    expect_identical(
        grep("^PG.Quantity_", names(d$protein_quant), value = TRUE),
        paste0("PG.Quantity_", d$metadata$Coding))
})

test_that("condition_order doubles as a filter", {
    d <- preprocess_diann(.diann_report(), condition_order = c("D", "A"),
                          verbose = FALSE)

    expect_identical(levels(d$metadata$R.Condition), c("D", "A"))
    expect_identical(unique(as.character(d$metadata$R.Condition)), c("D", "A"))
    expect_false(any(grepl("_B_", names(d$protein_quant), fixed = TRUE)))
    expect_identical(nrow(d$metadata), 8L)
})

test_that("sample_names reads the report as DIA-NN wrote it", {
    # The shipped example has already been renamed, so the two routes are
    # exercised against the same file by pretending the headers are opaque.
    d_named <- preprocess_diann(.diann_report(),
                                condition_order = c("A", "B", "D"),
                                verbose = FALSE)
    d_pos <- preprocess_diann(
        .diann_report(), condition_order = c("A", "B", "D"),
        sample_names = c(paste0("A_", 1:4), paste0("B_", 1:4),
                         paste0("D_", 1:4)),
        verbose = FALSE)

    expect_identical(d_pos$protein_quant, d_named$protein_quant)
    expect_identical(d_pos$protein_id,    d_named$protein_id)
    expect_identical(d_pos$metadata,      d_named$metadata)
})

test_that("an undeclared design fails with the columns listed in order", {
    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    df <- data.frame(Protein.Group = c("P1", "P2"),
                     Genes = c("G1", "G2"),
                     `run_one.raw` = c(1, 2),
                     `run_two.raw` = c(3, 4),
                     check.names = FALSE)
    utils::write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)

    expect_error(
        preprocess_diann(tmp, condition_order = "A", verbose = FALSE),
        "run_one.raw", fixed = TRUE)
    expect_error(
        preprocess_diann(tmp, condition_order = "A", verbose = FALSE),
        "sample_names")
})

test_that("input validation catches the usual mistakes", {
    f <- .diann_report()

    expect_error(preprocess_diann("no_such_file.tsv", condition_order = "A"),
                 "File not found")
    expect_error(preprocess_diann(f, condition_order = character(0)),
                 "non-empty character vector")
    expect_error(preprocess_diann(f, condition_order = "Z", verbose = FALSE),
                 "No run matches condition_order")
    expect_error(
        preprocess_diann(f, condition_order = "A",
                         sample_names = c("A_1", "A_2"), verbose = FALSE),
        "12 run columns")
    expect_error(
        preprocess_diann(f, condition_order = "A",
                         sample_names = rep("A_1", 12), verbose = FALSE),
        "duplicated")
    expect_error(
        preprocess_diann(f, condition_order = "A",
                         sample_names = c(paste0("A_", 1:11), "no_replicate"),
                         verbose = FALSE),
        "Could not parse these sample_names")

    # A file that is not a protein-group matrix at all
    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    utils::write.table(data.frame(a = 1, b = 2), tmp, sep = "\t",
                       row.names = FALSE, quote = FALSE)
    expect_error(preprocess_diann(tmp, condition_order = "A", verbose = FALSE),
                 "Protein.Group")
})

test_that("the metric families are exactly the ones DIA-NN can justify", {
    d <- preprocess_diann(.diann_report(), condition_order = c("A", "B", "D"),
                          verbose = FALSE)
    smp <- d$metadata$Coding

    # Fabricated all-NA, but present: both reactable loaders detect their
    # sample columns by grepping for these and abort without them.
    for (fam in c("PG.NrOfPrecursorsIdentified")) {
        cols <- paste0(fam, "_", smp)
        expect_true(all(cols %in% names(d$protein_id)))
        expect_true(all(vapply(d$protein_id[cols], is.numeric, logical(1))))
        expect_true(all(vapply(d$protein_id[cols], function(x) all(is.na(x)),
                               logical(1))))
    }
    cols <- paste0("PG.NrOfPrecursorsUsedForQuantification_", smp)
    expect_true(all(cols %in% names(d$protein_quant)))
    expect_true(all(vapply(d$protein_quant[cols], is.numeric, logical(1))))

    # Real data, from N.Sequences and N.Proteotypic.Sequences respectively.
    expect_false(anyNA(d$protein_id[[paste0("PG.NrOfStrippedSequencesIdentified_",
                                            smp[1])]]))
    expect_false(anyNA(
        d$protein_quant[[paste0("PG.NrOfStrippedSequencesUsedForQuantification_",
                                smp[1])]]))

    # Omitted outright: nothing downstream requires them, and an all-NA column
    # would draw a filter box that empties the table.
    expect_false(any(grepl("^PG.Coverage", names(d$protein_id))))
    expect_false(any(grepl("^PG.Cscore",   names(d$protein_id))))
    expect_false(any(grepl("^PG.Coverage", names(d$protein_quant))))
    expect_false(any(grepl("^PG.Cscore",   names(d$protein_quant))))

    # NA_real_, not NA: a logical column would come back from a .nadia file as
    # logical and break the identical() of the round trip.
    expect_type(d$protein_id$PG.MolecularWeight, "double")

    # protein_id and protein_quant share a row order, which write_nadia() needs.
    expect_identical(d$protein_id$PG.ProteinGroups,
                     d$protein_quant$PG.ProteinGroups)
})

test_that("metadata carries the one per-run count a matrix can give", {
    d <- preprocess_diann(.diann_report(), condition_order = c("A", "B", "D"),
                          verbose = FALSE)

    q <- as.matrix(d$protein_quant[paste0("PG.Quantity_", d$metadata$Coding)])
    expect_identical(d$metadata$R.ProteinGroupsIdentified,
                     as.numeric(colSums(!is.na(q))))
    # The D condition is the sparse one in this experiment.
    expect_lt(mean(d$metadata$R.ProteinGroupsIdentified[9:12]),
              mean(d$metadata$R.ProteinGroupsIdentified[1:8]))
})

test_that("the DIA-NN pipeline runs to completion", {
    d <- preprocess_diann(.diann_report(), condition_order = c("A", "B", "D"),
                          verbose = FALSE)
    res <- process_proteomics(d, verbose = FALSE)

    expect_s3_class(res, "proteomics_result")
    expect_length(res$comparisons, 3L)
    expect_gt(nrow(res$DEPs_results), 0)
    expect_false(anyNA(res$DEPs_results$logFC))

    # UniqPepts comes from N.Proteotypic.Sequences, so it is real and its
    # legitimate zeros must not have been turned into -Inf.
    up <- SummarizedExperiment::rowData(res$se_proc)$UniqPepts
    expect_true(all(is.finite(up)))
    expect_gt(max(up), 0)
})

test_that("the reactable tables open on a DIA-NN object", {
    # This is what the fabricated precursor families are for. The two loaders
    # abort without them, and no other format in the package exercises a metric
    # family that is entirely NA.
    d <- preprocess_diann(.diann_report(), condition_order = c("A", "B", "D"),
                          verbose = FALSE)

    expect_s3_class(protein_list_reactable(d$protein_id,
                                           metadata = d$metadata), "reactable")
    expect_s3_class(summary_list_widget(d$metadata), "shiny.tag.list")
})

test_that("preprocess_diann writes nothing without export_dir", {
    before <- list.files(tempdir(), recursive = TRUE)
    preprocess_diann(.diann_report(), condition_order = c("A", "B", "D"),
                     verbose = FALSE)
    expect_identical(list.files(tempdir(), recursive = TRUE), before)
})

test_that("the provenance attributes travel with the object", {
    d <- preprocess_diann(.diann_report(), condition_order = c("A", "B", "D"),
                          verbose = FALSE)

    expect_identical(names(d), c("metadata", "protein_id", "protein_quant"))
    expect_false(is.null(attr(d, "nadia_call")))
    src <- attr(d, "nadia_source")
    expect_identical(src$role, "report")
    expect_false(is.na(src$md5))
})
