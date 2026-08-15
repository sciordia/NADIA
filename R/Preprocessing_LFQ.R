# =============================================================================
# Proteome Discoverer LFQ Data Preprocessing
# =============================================================================
#
# Converts Proteome Discoverer LFQ protein exports into the same structure as
# `preprocess_spectronaut()` / `preprocess_tmt()`, so the downstream pipeline
# (Processing.R and the associated modules) consumes them unchanged.
#
# Unlike TMT (global per-protein metrics), an LFQ experiment carries PER-SAMPLE
# metrics in the same file, just like the Spectronaut report:
#   - `Abundance: <sample>`                       (intensity)
#   - `# PSMs (by Search Engine): <sample>`       (PSMs per sample)
#   - `# Peptides (by Search Engine): <sample>`   (peptides per sample)
#   - `Score Mascot: <sample>`                    (score per sample)
#
# The sample->condition mapping comes from the `Abundance:` suffixes, as in TMT
# and DIA-NN, or from a `Column`/`Condition` sheet when the columns are not
# named that way.
#
# All four families are mapped BY SAMPLE NAME, never by position, and that is
# not a stylistic choice: Proteome Discoverer writes them in different orders.
# Measured on a real export, `Score Mascot:` came out WT_1, MUT_1, WT_2, ...
# while `Abundance:` came out MUT_1, MUT_2, MUT_3, WT_1, ... Anything positional
# would have silently paired the score of one run with the intensity of another.
#
# Copyright 2025 Sergio Ciordia
# License: GPL-3
# =============================================================================

# --- Dependencies ---

# =============================================================================
# Internal Helper Functions
# =============================================================================

# .parse_gene_from_description() and .validate_pd_columns() live in R/utils.R:
# they used to be duplicated here and in Preprocessing_TMT.R.

#' Extract the species (OS=...) from the Description column
#' @description Takes the text between `OS=` and the next `XX=` token (e.g. OX=).
#' @noRd
.parse_species_from_description <- function(x) {
  m <- stringr::str_match(x, "OS=(.+?)\\s+[A-Za-z]+=")
  m[, 2]
}

#' Lookup that tolerates the usual variants of PD column names
#' @noRd
.pd_resolve_col <- function(df, ...) {
  candidates <- unlist(list(...), use.names = FALSE)
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

#' Extract a numeric vector from `df[[col]]`, or NA if the column is absent
#' @noRd
.pd_numeric_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) return(rep(NA_real_, nrow(df)))
  suppressWarnings(as.numeric(df[[col]]))
}

#' Extract a character vector from `df[[col]]`, or NA if the column is absent
#' @noRd
.pd_char_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) return(rep(NA_character_, nrow(df)))
  as.character(df[[col]])
}

#' Map the columns of a family (by prefix) to sample names
#' @description Runs `grep(prefix_regex, ...)`, strips the prefix to obtain the
#'   sample name and returns a named vector (sample name -> column name)
#'   restricted to the samples in the design. The mapping is BY NAME, so it is
#'   robust to ordering differences between column families.
#' @return Character vector named by sample (a subset of `samples`).
#' @noRd
.lfq_sample_cols <- function(df, prefix_regex, samples) {
  cols <- grep(prefix_regex, names(df), value = TRUE)
  if (length(cols) == 0) return(setNames(character(0), character(0)))
  smp <- sub(prefix_regex, "", cols)
  keep <- smp %in% samples
  setNames(cols[keep], smp[keep])
}

# =============================================================================
# Main Function
# =============================================================================

#' Preprocess Proteome Discoverer LFQ exports
#'
#' @description
#' Converts a Proteome Discoverer LFQ protein export (wide TSV format) into the
#' same output structure as `preprocess_spectronaut()` / `preprocess_tmt()`:
#' three data.frames (`metadata`, `protein_id`, `protein_quant`) ready for the
#' downstream pipeline (`process_proteomics()`).
#'
#' The experimental design comes from the `Abundance: <Condition>_<Replicate>`
#' column suffixes, exactly as in [preprocess_tmt()] and [preprocess_diann()].
#' For a report whose columns are not named that way, `annot_path` takes a sample
#' sheet instead.
#'
#' The intensities (`Abundance: <sample>`) and the three per-sample metric
#' families (`Score Mascot`, `# PSMs (by Search Engine)`,
#' `# Peptides (by Search Engine)`) are mapped BY SAMPLE NAME, never by position:
#' Proteome Discoverer writes the four families in different orders.
#'
#' @param file_path Path to the data TSV exported from Proteome Discoverer.
#' @param annot_path Path to a sample sheet with columns `Column` and
#'   `Condition`, where `Column` holds the text after `Abundance: `. If `NULL`
#'   (default), the design comes from the suffixes. Any other column in the sheet
#'   is ignored -- batch structure belongs in `covariate_df` of
#'   [process_proteomics()], see `vignette("batch-correction")`.
#' @param condition_order Character vector with the order of the conditions. If
#'   `NULL` (default), it is derived from the data (order of appearance). If
#'   supplied, it fixes the factor levels and discards samples whose condition is
#'   not in the list; conditions listed but absent are dropped too, with a
#'   warning, so that no empty factor level reaches the metadata.
#' @param export_dir Directory to export the TSVs to. `NULL` (default) = no
#'   export.
#' @param timestamp_suffix Logical. If `TRUE` (default), appends a timestamp to
#'   the exported files.
#' @param verbose Logical. If `TRUE` (default), shows progress messages.
#'
#' @return A list with class `c("lfq_data", "proteomics_data", "list")`
#'   containing `metadata`, `protein_id` and `protein_quant`.
#'
#'   The call that produced the object and a fingerprint of the files it was read
#'   from -- the report, and the sample sheet when there is one -- (path, size,
#'   modification time and MD5) travel with it as the attributes `nadia_call` and
#'   `nadia_source`. They
#'   are attributes rather than list elements so that the three-element structure
#'   above is unchanged; [write_nadia()] records them as the provenance of an
#'   analysis.
#'
#' @examples
#' # A trimmed Proteome Discoverer LFQ report ships with the package. Its columns
#' # follow the "Abundance: <Condition>_<Replicate>" convention, so the design
#' # needs no further argument.
#' report <- system.file("extdata", "nadia_lfq_report.tsv.gz", package = "NADIA")
#'
#' lfq <- preprocess_lfq(report, verbose = FALSE)
#' lfq
#' lfq$metadata[, c("Coding", "R.Condition", "R.Replicate")]
#'
#' # The same design, declared in a sheet instead. Use this route when the
#' # columns are not named by the convention; here it is the same experiment, so
#' # it produces the same sample table.
#' annot <- system.file("extdata", "nadia_lfq_annotation.tsv", package = "NADIA")
#' from_sheet <- preprocess_lfq(report, annot_path = annot, verbose = FALSE)
#' identical(from_sheet$metadata, lfq$metadata)
#'
#' @export
preprocess_lfq <- function(
    file_path,
    annot_path = NULL,
    condition_order = NULL,
    export_dir = NULL,
    timestamp_suffix = TRUE,
    verbose = TRUE
) {

  # --- Argument validation ---
  if (!file.exists(file_path)) stop("Data file not found: ", file_path)
  if (!is.null(condition_order) &&
      (length(condition_order) == 0 || !is.character(condition_order))) {
    stop("condition_order must be NULL or a non-empty character vector.")
  }

  # --- Reading (check.names = FALSE to preserve the PD headers) ---
  if (verbose) message("Reading data file: ", basename(file_path))
  df <- read.delim(file_path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
  .validate_pd_columns(df)

  # --- The design: the Abundance suffixes, or the sample sheet ---
  if (verbose && !is.null(annot_path)) {
    message("Reading annotation file: ", basename(annot_path))
  }
  design <- .pd_resolve_design(df, annot_path, verbose = verbose)
  from_sheet <- !is.null(annot_path)

  # --- Resolve the condition order and filter the design ---
  if (is.null(condition_order)) {
    condition_order <- unique(design$R.Condition)
  } else {
    keep <- design$R.Condition %in% condition_order
    if (!any(keep)) {
      stop("No condition in the ",
           if (from_sheet) "annotation file" else "report",
           " matches condition_order = c(",
           paste0("'", condition_order, "'", collapse = ", "), ").\n",
           "Conditions detected: ",
           paste(unique(design$R.Condition), collapse = ", "))
    }
    design <- design[keep, , drop = FALSE]
    condition_order <- .drop_absent_conditions(
      condition_order, design$R.Condition,
      source = if (from_sheet) "annotation file" else "report")
  }

  # --- Sample table, sorted by condition and then by replicate ---
  design$R.Condition <- factor(design$R.Condition, levels = condition_order,
                               ordered = TRUE)
  design <- design[order(design$R.Condition, design$R.Replicate), , drop = FALSE]
  coding_levels <- design$Coding

  # ==========================================================================
  # metadata
  # ==========================================================================
  if (verbose) message("Generating sample metadata...")
  run_summary <- data.frame(
    R.FileName   = design$abundance_col,
    R.Condition  = design$R.Condition,
    R.Replicate  = design$R.Replicate,
    Coding       = coding_levels,
    stringsAsFactors = FALSE
  )
  rownames(run_summary) <- run_summary$Coding

  # ==========================================================================
  # Base protein information
  # ==========================================================================
  if (verbose) message("Processing the base protein information...")

  col_mw     <- .pd_resolve_col(df, "MW [kDa]", "MW (kDa)", "MW")
  col_pi     <- .pd_resolve_col(df, "calc. pI", "calc pI", "Calculated pI")
  col_master <- .pd_resolve_col(df, "Master")
  col_pgids  <- .pd_resolve_col(df, "Protein Group IDs")
  col_fdr    <- .pd_resolve_col(df, "Protein FDR Confidence: Combined",
                                    "Protein FDR Confidence")
  # Proteome Discoverer 3.3 renamed this field to "Exp. Protein q-value".
  col_qvalue <- .pd_resolve_col(df, "Exp. Protein q-value: Combined",
                                    "Exp. Protein q-value",
                                    "Exp. q-value: Combined", "Exp. q-value")
  col_pep    <- .pd_resolve_col(df, "Sum PEP Score")
  col_psms   <- .pd_resolve_col(df, "# PSMs", "Number of PSMs")
  col_pepts  <- .pd_resolve_col(df, "# Peptides", "Number of Peptides")
  col_uniq   <- .pd_resolve_col(df, "# Unique Peptides", "Number of Unique Peptides")
  col_cov    <- .pd_resolve_col(df, "Coverage [%]", "Coverage (%)", "Coverage")

  base_info <- data.frame(
    PG.ProteinGroups       = as.character(df$Accession),
    PG.ProteinDescriptions = as.character(df$Description),
    PG.Genes               = .parse_gene_from_description(df$Description),
    PG.MolecularWeight     = .pd_numeric_or_na(df, col_mw),
    stringsAsFactors = FALSE
  )

  lfq_extras <- data.frame(
    calc.pI                = .pd_numeric_or_na(df, col_pi),
    Master                 = .pd_char_or_na(df, col_master),
    Protein.Group.IDs      = .pd_char_or_na(df, col_pgids),
    Protein.FDR.Confidence = .pd_char_or_na(df, col_fdr),
    Exp.q.value            = .pd_numeric_or_na(df, col_qvalue),
    NrUniquePeptides       = .pd_numeric_or_na(df, col_uniq),
    Species                = .parse_species_from_description(df$Description),
    stringsAsFactors = FALSE
  )

  # Global metrics
  g_psms <- .pd_numeric_or_na(df, col_psms)
  g_pept <- .pd_numeric_or_na(df, col_pepts)
  g_cov  <- .pd_numeric_or_na(df, col_cov)
  g_pep  <- .pd_numeric_or_na(df, col_pep)

  # --- Helper: per-sample matrix (mapped by name) for one column family ---
  sample_matrix <- function(prefix_regex, out_prefix) {
    map <- .lfq_sample_cols(df, prefix_regex, coding_levels)
    mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(coding_levels),
                  dimnames = list(NULL, paste0(out_prefix, "_", coding_levels)))
    for (i in seq_along(coding_levels)) {
      col <- map[[coding_levels[i]]]
      if (!is.null(col) && !is.na(col)) {
        mat[, i] <- suppressWarnings(as.numeric(df[[col]]))
      }
    }
    as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  }
  # Replicate a global metric across samples (to get the wide shape)
  wide_global <- function(values, out_prefix) {
    mat <- matrix(rep(values, length(coding_levels)),
                  nrow = length(values), ncol = length(coding_levels))
    colnames(mat) <- paste0(out_prefix, "_", coding_levels)
    as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # Per-sample column families
  psms_by <- sample_matrix("^# PSMs \\(by Search Engine\\):\\s*",
                           "PG.NrOfPrecursorsUsedForQuantification")
  pept_by <- sample_matrix("^# Peptides \\(by Search Engine\\):\\s*",
                           "PG.NrOfStrippedSequencesUsedForQuantification")
  mascot_by <- sample_matrix("^Score Mascot:\\s*", "PG.Cscore.RunWise")

  abund_mat <- df[, design$abundance_col, drop = FALSE]
  abund_mat[] <- lapply(abund_mat, function(v) suppressWarnings(as.numeric(v)))
  names(abund_mat) <- paste0("PG.Quantity_", coding_levels)

  # ==========================================================================
  # protein_ID  (per-sample identification metrics wherever they exist)
  # ==========================================================================
  if (verbose) message("Processing protein_ID...")

  id_psms <- psms_by; names(id_psms) <- paste0("PG.NrOfPrecursorsIdentified_", coding_levels)
  id_pept <- pept_by; names(id_pept) <- paste0("PG.NrOfStrippedSequencesIdentified_", coding_levels)
  id_cov  <- wide_global(g_cov, "PG.Coverage")
  id_score <- mascot_by  # PG.Cscore.RunWise_<coding>

  protein_ID <- cbind(base_info, lfq_extras, id_psms, id_pept, id_cov, id_score)

  static_cols  <- names(base_info)
  extras_cols  <- names(lfq_extras)
  metric_order <- c("PG.NrOfPrecursorsIdentified",
                    "PG.NrOfStrippedSequencesIdentified",
                    "PG.Coverage", "PG.Cscore.RunWise")
  runwise_order <- unlist(lapply(metric_order,
                                 function(m) paste0(m, "_", coding_levels)))
  protein_ID <- protein_ID[, c(static_cols, extras_cols,
                               intersect(runwise_order, names(protein_ID))),
                           drop = FALSE]
  protein_ID <- protein_ID[order(protein_ID$PG.ProteinGroups), , drop = FALSE]
  rownames(protein_ID) <- NULL
  if (anyDuplicated(protein_ID$PG.ProteinGroups) > 0) {
    stop("Integrity check failed: protein_ID contains duplicated Accession values.")
  }

  # ==========================================================================
  # protein_QUANT
  # ==========================================================================
  if (verbose) message("Processing protein_QUANT...")

  global_metrics <- data.frame(
    PG.NrOfPrecursorsIdentified.Global        = g_psms,
    PG.NrOfStrippedSequencesIdentified.Global = g_pept,
    PG.Coverage.Global                        = g_cov,
    PG.Cscore                                 = g_pep,
    stringsAsFactors = FALSE
  )

  protein_QUANT <- cbind(base_info, global_metrics, psms_by, pept_by, abund_mat)

  static_cols2 <- names(base_info)
  global_cols  <- names(global_metrics)
  q_metric_order <- c("PG.NrOfPrecursorsUsedForQuantification",
                      "PG.NrOfStrippedSequencesUsedForQuantification",
                      "PG.Quantity")
  q_runwise <- unlist(lapply(q_metric_order,
                             function(m) paste0(m, "_", coding_levels)))
  protein_QUANT <- protein_QUANT[, c(static_cols2, global_cols,
                                     intersect(q_runwise, names(protein_QUANT))),
                                 drop = FALSE]
  protein_QUANT <- protein_QUANT[order(protein_QUANT$PG.ProteinGroups), , drop = FALSE]
  rownames(protein_QUANT) <- NULL
  if (anyDuplicated(protein_QUANT$PG.ProteinGroups) > 0) {
    stop("Integrity check failed: protein_QUANT contains duplicated Accession values.")
  }

  # ==========================================================================
  # Optional export
  # ==========================================================================
  if (!is.null(export_dir)) {
    if (verbose) message("Exporting files to: ", export_dir)
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
    suffix <- if (timestamp_suffix) paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S")) else ""
    readr::write_tsv(run_summary,
                     file.path(export_dir, paste0("Metadata", suffix, ".tsv")), na = "")
    readr::write_tsv(protein_ID,
                     file.path(export_dir, paste0("Protein_ID", suffix, ".tsv")), na = "")
    readr::write_tsv(protein_QUANT,
                     file.path(export_dir, paste0("Protein_QUANT", suffix, ".tsv")), na = "")
    if (verbose) message("Files exported successfully.")
  }

  # ==========================================================================
  # Result
  # ==========================================================================
  if (verbose) {
    message("Processing complete:\n",
            "  - Samples: ", nrow(run_summary), "\n",
            "  - Proteins (ID): ", nrow(protein_ID), "\n",
            "  - Proteins (QUANT): ", nrow(protein_QUANT))
  }

  result <- list(
    metadata      = as.data.frame(run_summary),
    protein_id    = as.data.frame(protein_ID),
    protein_quant = as.data.frame(protein_QUANT)
  )
  class(result) <- c("lfq_data", "proteomics_data", "list")
  result <- .nadia_stamp(result, match.call(),
                         c(report = file_path, annotation = annot_path))
  return(result)
}

# =============================================================================
# Methods for the lfq_data class
# =============================================================================

#' Print a summary of an LFQ preprocessing result
#'
#' Reports the number of samples, protein groups and conditions, so that the
#' shape of the experiment can be checked at a glance before processing it.
#'
#' @param x An `lfq_data` object, as returned by [preprocess_lfq()].
#' @param ... Ignored, present for compatibility with the `print` generic.
#' @return `x`, invisibly. Called for the summary it prints.
#'
#' @examples
#' lfq <- preprocess_lfq(
#'   system.file("extdata", "nadia_lfq_report.tsv.gz", package = "NADIA"),
#'   annot_path = system.file("extdata", "nadia_lfq_annotation.tsv",
#'                            package = "NADIA"),
#'   verbose = FALSE)
#' print(lfq)
#'
#' @export
print.lfq_data <- function(x, ...) {
  cat("Preprocessed LFQ (Proteome Discoverer) data\n")
  cat("------------------------------------------\n")
  cat("Samples (metadata):", nrow(x$metadata), "\n")
  cat("Proteins (ID):", nrow(x$protein_id), "\n")
  cat("Proteins (QUANT):", nrow(x$protein_quant), "\n")
  cat("\nConditions:",
      paste(levels(x$metadata$R.Condition) %||%
              unique(x$metadata$R.Condition), collapse = ", "), "\n")
  invisible(x)
}
