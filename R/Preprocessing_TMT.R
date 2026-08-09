# =============================================================================
# Proteome Discoverer TMT Data Preprocessing
# =============================================================================
#
# Converts Proteome Discoverer protein exports (TMT/TMTpro) into the same
# structure as `preprocess_spectronaut()`, so the downstream pipeline
# (Processing.R and the associated modules) consumes them unchanged.
#
# Copyright 2025 Sergio Ciordia
# License: GPL-3
# =============================================================================

# --- Dependencies ---

# =============================================================================
# Internal Helper Functions
# =============================================================================

# .parse_gene_from_description() and .validate_pd_columns() live in R/utils.R:
# they used to be duplicated here and in Preprocessing_LFQ.R.
#
# .parse_abundance_columns() also lives in R/utils.R: preprocess_diann() reads
# the same `Abundance: <Condition>_<Replicate>` convention, so the parser is
# shared rather than copied a third time.

#' Lookup that tolerates the usual variants of PD column names
#' @description PD exports columns with special characters (spaces, #, %,
#'   brackets, ":"). This function tries several variants and returns the first
#'   one that exists, or `NA_character_` if none matches.
#' @noRd
.tmt_resolve_col <- function(df, ...) {
  candidates <- unlist(list(...), use.names = FALSE)
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

#' Extract a numeric vector from `df[[col]]`, or NA if the column is absent
#' @noRd
.tmt_numeric_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[col]]))
}

#' Extract a character vector from `df[[col]]`, or NA if the column is absent
#' @noRd
.tmt_char_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) {
    return(rep(NA_character_, nrow(df)))
  }
  as.character(df[[col]])
}

# =============================================================================
# Main Function
# =============================================================================

#' Preprocess Proteome Discoverer TMT exports
#'
#' @description
#' Converts a Proteome Discoverer protein export (wide TSV format, with
#' `Abundance: <Condition>_<Replicate>` columns) into the same output structure
#' as `preprocess_spectronaut()`: three data.frames (`metadata`, `protein_id`,
#' `protein_quant`) ready for the downstream pipeline (`process_proteomics()`
#' and the associated modules).
#'
#' The global PD identification metrics (`# PSMs`, `# Peptides`,
#' `# Unique Peptides`, `Coverage [%]`) are replicated into per-channel columns
#' to fit the wide Spectronaut contract.
#'
#' @param file_path Path to the TSV exported from Proteome Discoverer.
#' @param condition_order Character vector with the order of the experimental
#'   conditions (e.g. `c("A","B","C","D","IS")`). Only the channels whose
#'   condition is in this vector are kept -- handy for excluding Internal
#'   Standards by omitting `"IS"`.
#' @param export_dir Directory to export the TSV files to. If `NULL` (default),
#'   no files are exported.
#' @param timestamp_suffix Logical. If `TRUE` (default), appends a timestamp to
#'   the names of the exported files.
#' @param verbose Logical. If `TRUE` (default), shows progress messages.
#'
#' @return A list with class `c("tmt_data", "proteomics_data", "list")`
#'   containing:
#'   \describe{
#'     \item{metadata}{Data frame with one record per channel/sample}
#'     \item{protein_id}{Data frame with the identification metrics per protein
#'       (wide format, global metrics replicated per channel)}
#'     \item{protein_quant}{Data frame with the quantification metrics per
#'       protein (includes `PG.Quantity_<Coding>`)}
#'   }
#'
#'   The call that produced the object and a fingerprint of the file it was read
#'   from (path, size, modification time and MD5) travel with it as the
#'   attributes `nadia_call` and `nadia_source`. They are attributes rather than
#'   list elements so that the three-element structure above is unchanged;
#'   [write_nadia()] records them as the provenance of an analysis.
#'
#' @examples
#' # A trimmed Proteome Discoverer TMTpro report ships with the package
#' report <- system.file("extdata", "nadia_tmt_report.tsv.gz", package = "NADIA")
#'
#' # condition_order is required, and doubles as a filter: the report also holds
#' # an "IS" channel (internal standards) that is left out by not listing it
#' tmt <- preprocess_tmt(report,
#'                       condition_order = c("A", "B", "C", "D"),
#'                       verbose = FALSE)
#' tmt
#' table(tmt$metadata$R.Condition)
#'
#' @export
preprocess_tmt <- function(
    file_path,
    condition_order,
    export_dir = NULL,
    timestamp_suffix = TRUE,
    verbose = TRUE
) {

  # --- Argument validation ---
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  if (length(condition_order) == 0 || !is.character(condition_order)) {
    stop("condition_order must be a non-empty character vector.")
  }

  # --- Reading the file (check.names = FALSE to preserve the PD headers) ---
  if (verbose) message("Reading file: ", basename(file_path))

  df <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # --- Validate the minimum set of columns ---
  .validate_pd_columns(df)

  # --- Parse the Abundance columns ---
  abund <- .parse_abundance_columns(df)

  # Filter the channels by condition_order (drops e.g. "IS" when not listed)
  keep <- abund$R.Condition %in% condition_order
  if (!any(keep)) {
    stop(
      "No Abundance channel matches condition_order = c(",
      paste0("'", condition_order, "'", collapse = ", "), ").\n",
      "Conditions detected in the file: ",
      paste(unique(abund$R.Condition), collapse = ", ")
    )
  }
  abund <- abund[keep, , drop = FALSE]

  # Final ordering by condition_order and replicate
  abund$R.Condition <- factor(abund$R.Condition,
                              levels = condition_order, ordered = TRUE)
  abund <- abund[order(abund$R.Condition, abund$R.Replicate), , drop = FALSE]
  coding_levels <- abund$Coding

  # --- Resolve the optional PD columns (tolerant name matching) ---
  col_mw       <- .tmt_resolve_col(df, "MW [kDa]", "MW (kDa)", "MW")
  col_pi       <- .tmt_resolve_col(df, "calc. pI", "calc pI", "Calculated pI")
  col_master   <- .tmt_resolve_col(df, "Master")
  col_pgids    <- .tmt_resolve_col(df, "Protein Group IDs")
  col_fdr      <- .tmt_resolve_col(df, "Protein FDR Confidence: Combined",
                                       "Protein FDR Confidence")
  col_qvalue   <- .tmt_resolve_col(df, "Exp. q-value: Combined",
                                       "Exp. q-value")
  col_pep      <- .tmt_resolve_col(df, "Sum PEP Score")
  col_mascot   <- .tmt_resolve_col(df, "Score Mascot: Mascot", "Score: Mascot")
  col_sequest  <- .tmt_resolve_col(df, "Score Sequest HT: Sequest HT",
                                       "Score: Sequest HT")
  col_psms     <- .tmt_resolve_col(df, "# PSMs", "Number of PSMs")
  col_pepts    <- .tmt_resolve_col(df, "# Peptides", "Number of Peptides")
  col_uniq     <- .tmt_resolve_col(df, "# Unique Peptides",
                                       "Number of Unique Peptides")
  col_cov      <- .tmt_resolve_col(df, "Coverage [%]", "Coverage (%)",
                                       "Coverage")

  # ==========================================================================
  # metadata
  # ==========================================================================
  if (verbose) message("Generating channel metadata...")

  run_summary <- data.frame(
    R.FileName  = abund$abundance_col,
    R.Condition = abund$R.Condition,
    R.Replicate = abund$R.Replicate,
    Coding      = abund$Coding,
    stringsAsFactors = FALSE
  )
  rownames(run_summary) <- run_summary$Coding

  # ==========================================================================
  # Base protein information (shared between protein_ID and protein_QUANT)
  # ==========================================================================
  if (verbose) message("Processing the base protein information...")

  base_info <- data.frame(
    PG.ProteinGroups       = as.character(df$Accession),
    PG.ProteinDescriptions = as.character(df$Description),
    PG.Genes               = .parse_gene_from_description(df$Description),
    PG.MolecularWeight     = .tmt_numeric_or_na(df, col_mw),
    stringsAsFactors = FALSE
  )

  # TMT-specific extras (kept after MW so as not to break downstream code)
  tmt_extras <- data.frame(
    calc.pI                = .tmt_numeric_or_na(df, col_pi),
    Master                 = .tmt_char_or_na(df, col_master),
    Protein.Group.IDs      = .tmt_char_or_na(df, col_pgids),
    Protein.FDR.Confidence = .tmt_char_or_na(df, col_fdr),
    Exp.q.value            = .tmt_numeric_or_na(df, col_qvalue),
    Score.Mascot           = .tmt_numeric_or_na(df, col_mascot),
    Score.SequestHT        = .tmt_numeric_or_na(df, col_sequest),
    NrUniquePeptides       = .tmt_numeric_or_na(df, col_uniq),
    stringsAsFactors = FALSE
  )

  # Global metrics that will be replicated per channel
  g_psms <- .tmt_numeric_or_na(df, col_psms)
  g_pept <- .tmt_numeric_or_na(df, col_pepts)
  g_cov  <- .tmt_numeric_or_na(df, col_cov)
  g_pep  <- .tmt_numeric_or_na(df, col_pep)
  g_uniq <- .tmt_numeric_or_na(df, col_uniq)

  # ==========================================================================
  # protein_ID  (global metrics replicated per channel to get the wide shape)
  # ==========================================================================
  if (verbose) message("Processing protein_ID...")

  mk_wide <- function(values, metric_name) {
    mat <- matrix(rep(values, length(coding_levels)),
                  nrow = length(values), ncol = length(coding_levels))
    colnames(mat) <- paste0(metric_name, "_", coding_levels)
    as.data.frame(mat, stringsAsFactors = FALSE)
  }

  id_psms <- mk_wide(g_psms, "PG.NrOfPrecursorsIdentified")
  id_pept <- mk_wide(g_pept, "PG.NrOfStrippedSequencesIdentified")
  id_cov  <- mk_wide(g_cov,  "PG.Coverage")
  id_pep  <- mk_wide(g_pep,  "PG.Cscore.RunWise")

  protein_ID <- cbind(base_info, tmt_extras, id_psms, id_pept, id_cov, id_pep)

  # Order the columns: static + extras + metrics in a stable order
  static_cols  <- c("PG.ProteinGroups", "PG.ProteinDescriptions",
                    "PG.Genes", "PG.MolecularWeight")
  extras_cols  <- names(tmt_extras)
  metric_order <- c("PG.NrOfPrecursorsIdentified",
                    "PG.NrOfStrippedSequencesIdentified",
                    "PG.Coverage", "PG.Cscore.RunWise")
  runwise_order <- unlist(lapply(metric_order,
                                 function(m) paste0(m, "_", coding_levels)))
  final_cols <- c(static_cols, extras_cols,
                  intersect(runwise_order, names(protein_ID)))
  protein_ID <- protein_ID[, final_cols, drop = FALSE]
  protein_ID <- protein_ID[order(protein_ID$PG.ProteinGroups), , drop = FALSE]
  rownames(protein_ID) <- NULL

  # Check uniqueness
  if (anyDuplicated(protein_ID$PG.ProteinGroups) > 0) {
    stop("Integrity check failed: protein_ID contains duplicated Accession values. ",
         "Check the input data.")
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

  # Per-channel metrics: NrOfPrecursorsUsedForQuantification (replicated PSMs),
  # NrOfStrippedSequencesUsedForQuantification (replicated Unique Peptides),
  # Quantity (Abundance: <Coding>)
  q_psms <- mk_wide(g_psms, "PG.NrOfPrecursorsUsedForQuantification")
  q_uniq <- mk_wide(g_uniq, "PG.NrOfStrippedSequencesUsedForQuantification")

  abund_mat <- df[, abund$abundance_col, drop = FALSE]
  abund_mat[] <- lapply(abund_mat, function(v) suppressWarnings(as.numeric(v)))
  names(abund_mat) <- paste0("PG.Quantity_", coding_levels)

  protein_QUANT <- cbind(base_info, global_metrics, q_psms, q_uniq, abund_mat)

  static_cols2 <- c("PG.ProteinGroups", "PG.ProteinDescriptions",
                    "PG.Genes", "PG.MolecularWeight")
  global_cols  <- c("PG.NrOfPrecursorsIdentified.Global",
                    "PG.NrOfStrippedSequencesIdentified.Global",
                    "PG.Coverage.Global", "PG.Cscore")
  q_metric_order <- c("PG.NrOfPrecursorsUsedForQuantification",
                      "PG.NrOfStrippedSequencesUsedForQuantification",
                      "PG.Quantity")
  q_runwise <- unlist(lapply(q_metric_order,
                             function(m) paste0(m, "_", coding_levels)))
  final_cols2 <- c(static_cols2, global_cols,
                   intersect(q_runwise, names(protein_QUANT)))
  protein_QUANT <- protein_QUANT[, final_cols2, drop = FALSE]
  protein_QUANT <- protein_QUANT[order(protein_QUANT$PG.ProteinGroups), ,
                                 drop = FALSE]
  rownames(protein_QUANT) <- NULL

  if (anyDuplicated(protein_QUANT$PG.ProteinGroups) > 0) {
    stop("Integrity check failed: protein_QUANT contains duplicated Accession values. ",
         "Check the input data.")
  }

  # ==========================================================================
  # Optional export
  # ==========================================================================
  if (!is.null(export_dir)) {
    if (verbose) message("Exporting files to: ", export_dir)

    if (!dir.exists(export_dir)) {
      dir.create(export_dir, recursive = TRUE)
    }

    suffix <- if (timestamp_suffix) {
      paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    } else {
      ""
    }

    readr::write_tsv(
      run_summary,
      file.path(export_dir, paste0("Metadata", suffix, ".tsv")),
      na = ""
    )
    readr::write_tsv(
      protein_ID,
      file.path(export_dir, paste0("Protein_ID", suffix, ".tsv")),
      na = ""
    )
    readr::write_tsv(
      protein_QUANT,
      file.path(export_dir, paste0("Protein_QUANT", suffix, ".tsv")),
      na = ""
    )

    if (verbose) message("Files exported successfully.")
  }

  # ==========================================================================
  # Build the result
  # ==========================================================================
  if (verbose) {
    message(
      "Processing complete:\n",
      "  - Channels: ", nrow(run_summary), "\n",
      "  - Proteins (ID): ", nrow(protein_ID), "\n",
      "  - Proteins (QUANT): ", nrow(protein_QUANT)
    )
  }

  result <- list(
    metadata      = as.data.frame(run_summary),
    protein_id    = as.data.frame(protein_ID),
    protein_quant = as.data.frame(protein_QUANT)
  )

  class(result) <- c("tmt_data", "proteomics_data", "list")
  result <- .nadia_stamp(result, match.call(), c(report = file_path))
  return(result)
}

# =============================================================================
# Methods for the tmt_data class
# =============================================================================

#' Print a summary of a TMT preprocessing result
#'
#' Reports the number of channels, protein groups and conditions, so that the
#' shape of the experiment can be checked at a glance before processing it.
#'
#' @param x A `tmt_data` object, as returned by [preprocess_tmt()].
#' @param ... Ignored, present for compatibility with the `print` generic.
#' @return `x`, invisibly. Called for the summary it prints.
#'
#' @examples
#' tmt <- preprocess_tmt(
#'   system.file("extdata", "nadia_tmt_report.tsv.gz", package = "NADIA"),
#'   condition_order = c("A", "B", "C", "D"),
#'   verbose = FALSE)
#' print(tmt)
#'
#' @export
print.tmt_data <- function(x, ...) {
  cat("Preprocessed TMT (Proteome Discoverer) data\n")
  cat("------------------------------------------\n")
  cat("Channels (metadata):", nrow(x$metadata), "\n")
  cat("Proteins (ID):", nrow(x$protein_id), "\n")
  cat("Proteins (QUANT):", nrow(x$protein_quant), "\n")
  cat("\nConditions:",
      paste(levels(x$metadata$R.Condition) %||%
              unique(x$metadata$R.Condition), collapse = ", "),
      "\n")
  invisible(x)
}
