# =============================================================================
# NADIA: an end-to-end analysis, written to be adapted
# =============================================================================
#
# This is a template, not a turnkey script. Every block marked EDIT encodes a
# decision that belongs to the analyst. Read inst/skill/reference/ before
# changing them, and never report the defaults as if they were a conclusion.
#
# As shipped it runs on the bundled example dataset, so it works unmodified and
# you can see the shape of the output before pointing it at real data.
# =============================================================================

library(NADIA)

# --- EDIT 1: the input ------------------------------------------------------
# Replace with your own report. condition_order sets the factor levels, and so
# the direction of every contrast: the first level is the reference.
#
#   prep <- preprocess_spectronaut("report.tsv", condition_order = c("A", "B", "D"))
#   prep <- preprocess_diann("report.pg_matrix.tsv", condition_order = c("A", "B", "D"))
#   prep <- preprocess_tmt("tmt_report.txt", condition_order = c("A", "B", "D"))
#   prep <- preprocess_lfq("lfq_report.txt", annot_path = "samples.tsv")

data(nadia_dia)
prep <- nadia_dia

out_dir <- file.path(tempdir(), "nadia_run")   # EDIT: a real project directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Check the import before going further ----------------------------------
message("Samples per condition:")
print(table(prep$metadata$R.Condition))

quant <- prep$protein_quant[, grep("PG.Quantity", names(prep$protein_quant))]
if (ncol(quant)) {
    q <- as.matrix(quant)
    message(sprintf("Protein groups: %d | missing before imputation: %.1f %%",
                    nrow(prep$protein_quant),
                    100 * mean(is.na(q) | q == 0)))
}

# --- EDIT 2: the methods ----------------------------------------------------
# The defaults below are a starting point. Justify them with the metrics
# section further down, or with a spike-in benchmark, before reporting.
norm_method <- "cycloess"
imp_method  <- "combo"
mar_method  <- "Impseqrob"
mnar_method <- "min"
de_method   <- "limma"

alpha           <- 0.05
logFC_threshold <- 0

# --- Run the pipeline -------------------------------------------------------
res <- process_proteomics(
    prep,
    norm_method     = norm_method,
    imp_method      = imp_method,
    mar_method      = mar_method,
    mnar_method     = mnar_method,
    de_method       = de_method,
    alpha           = alpha,
    logFC_threshold = logFC_threshold,

    # EDIT 3: batch correction. Off unless you have a real batch structure.
    # Pass the biological variable in batch_covariates or ComBat can erase it.
    #   covariate_df     = my_batch_table,
    #   batch_correct    = TRUE,
    #   batch_column     = "Batch",
    #   batch_covariates = "Condition",

    # EDIT 4: contrasts. NULL = all pairwise; or comparisons = c("B-A"),
    # or control = "A" for everything against one condition.
    comparisons = NULL,

    verbose = FALSE
)

imputed_assay <- utils::tail(SummarizedExperiment::assayNames(res$se_proc), 1)
message(sprintf("Assays: %s",
                paste(SummarizedExperiment::assayNames(res$se_proc), collapse = ", ")))

# --- What came out ----------------------------------------------------------
message("\nClassification per comparison:")
print(table(res$DEPs_results$Comparison, res$DEPs_results$Change))

# How much of the significant set rests on imputed values? Report this.
sig <- res$DEPs_results[res$DEPs_results$Change != "No Change", ]
message("\nMissingness within the compared conditions (MissComp) of the significant set:")
print(summary(sig$MissComp))

# --- EDIT 5: justify the method choice --------------------------------------
# Without a ground truth, rank on data-derived metrics. Slow; enable when you
# are choosing methods rather than re-running a settled pipeline.
run_metrics <- FALSE
if (run_metrics) {
    nm <- normalization_metrics(res$se_proc,
                                output_dir = file.path(out_dir, "norm_metrics"),
                                verbose = FALSE)
    print(utils::head(nm$final_rank))

    im <- imputation_metrics(res$se_proc, assay_name = norm_method,
                             output_dir = file.path(out_dir, "imp_metrics"),
                             verbose = FALSE)
    print(utils::head(im$metrics_table))
}

# With a spike-in, measure recovery of the expected changes instead:
#   bm <- benchmarking_proteomics(res$DEPs_results, expected_values = expected,
#                                 species_df = species_map,
#                                 output_dir = file.path(out_dir, "benchmark"))

# --- EDIT 6: profiles across conditions (optional) --------------------------
run_pattern_profiler <- FALSE
if (run_pattern_profiler &&
    all(vapply(c("Mfuzz", "Biobase", "e1071"), requireNamespace, logical(1),
               quietly = TRUE))) {
    # assay_name is required in practice: its default names an assay that no
    # current pipeline produces.
    pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
                                    assay_name = imputed_assay,
                                    filter_mode = "any", seed = 123,
                                    verbose = FALSE)
    message(sprintf("Optimal number of clusters: %d", pp$optimal_c))

    # DE statistics joined to the dominant cluster: the input for enrichment.
    # A left join from DEPs_results: proteins without a cluster stay, as the
    # background.
    functional_input <- deps_with_clusters(res$DEPs_results, pp,
                                           assignment = "primary")
}

# --- EDIT 7: outputs --------------------------------------------------------
# Tables. export_format is "tsv", "parquet" or "both".
invisible(process_proteomics(prep,
                             norm_method = norm_method, imp_method = imp_method,
                             mar_method = mar_method, mnar_method = mnar_method,
                             de_method = de_method,
                             export_dir = file.path(out_dir, "tables"),
                             export_format = "tsv", verbose = FALSE))

# One archive holding the analysis, its parameters and its provenance.
# Never keep an open .nadia inside iCloud Drive, Dropbox or OneDrive.
if (requireNamespace("duckdb", quietly = TRUE)) {
    nadia_file <- file.path(out_dir, "analysis.nadia")
    write_nadia(nadia_file, preprocessing = prep, result = res,
                overwrite = TRUE, verbose = FALSE)
    message(sprintf("Archive written: %s (%.1f MB)",
                    nadia_file, file.size(nadia_file) / 1024^2))
}

# --- Record what produced the numbers ---------------------------------------
message("\nSettings used:")
print(unlist(res$parameters[c("norm_method", "imp_method", "mar_method",
                              "mnar_method", "de_method", "alpha",
                              "logFC_threshold")]))

message(sprintf("\nOutputs under: %s", out_dir))
