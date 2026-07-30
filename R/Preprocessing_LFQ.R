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
# The sample->condition mapping is taken from the annotation file (`_Annot`),
# with columns `Column`, `Condition` (mandatory) and `Experiment` (optional).
# Intensities are mapped BY SAMPLE NAME (not by position): the order of the
# `Abundance:` columns may differ from that of the metric columns.
#
# Copyright 2025 Sergio Ciordia
# Licensed under MIT
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

#' Validate the annotation file (experimental design)
#' @noRd
.validate_lfq_annot <- function(annot) {
  required <- c("Column", "Condition")
  missing <- setdiff(required, names(annot))
  if (length(missing) > 0) {
    stop(
      "Required columns missing from the annotation file (_Annot):\n  - ",
      paste(missing, collapse = "\n  - "),
      "\nThe _Annot must have at least 'Column' and 'Condition'."
    )
  }
  if (anyDuplicated(annot$Column) > 0) {
    stop("The annotation file has duplicated values in 'Column'.")
  }
  invisible(TRUE)
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
#' The experimental design (sample -> condition) is taken from the annotation
#' file `annot_path` (columns `Column`, `Condition`, and optionally
#' `Experiment`). The intensities (`Abundance: <sample>`) and the per-sample
#' metrics (`# PSMs (by Search Engine)`, `# Peptides (by Search Engine)`,
#' `Score Mascot`) are mapped BY SAMPLE NAME.
#'
#' @param file_path Path to the data TSV exported from Proteome Discoverer.
#' @param annot_path Path to the annotation TSV with columns `Column`,
#'   `Condition` (mandatory) and `Experiment` (optional).
#' @param condition_order Character vector with the order of the conditions. If
#'   `NULL` (default), it is derived from `Condition` (order of appearance). If
#'   supplied, it fixes the factor levels and discards samples whose condition is
#'   not in the list.
#' @param export_dir Directory to export the TSVs to. `NULL` (default) = no
#'   export.
#' @param timestamp_suffix Logical. If `TRUE` (default), appends a timestamp to
#'   the exported files.
#' @param verbose Logical. If `TRUE` (default), shows progress messages.
#'
#' @return A list with class `c("lfq_data", "proteomics_data", "list")`
#'   containing `metadata`, `protein_id` and `protein_quant`.
#'
#' @examples
#' \dontrun{
#' result <- preprocess_lfq(
#'   file_path  = "data-raw/20260710_AIturrate_2659_LFQ_QUANT_onlyRAW.tsv",
#'   annot_path = "data-raw/20260710_AIturrate_2659_LFQ_QUANT_onlyRAW_Annot.tsv",
#'   export_dir = "./results"
#' )
#' print(result)
#' }
#'
#' @export
preprocess_lfq <- function(
    file_path,
    annot_path,
    condition_order = NULL,
    export_dir = NULL,
    timestamp_suffix = TRUE,
    verbose = TRUE
) {

  # --- Argument validation ---
  if (!file.exists(file_path)) stop("Data file not found: ", file_path)
  if (!file.exists(annot_path)) stop("Annotation file not found: ", annot_path)
  if (!is.null(condition_order) &&
      (length(condition_order) == 0 || !is.character(condition_order))) {
    stop("condition_order must be NULL or a non-empty character vector.")
  }

  # --- Reading (check.names = FALSE to preserve the PD headers) ---
  if (verbose) message("Reading data file: ", basename(file_path))
  df <- read.delim(file_path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
  .validate_pd_columns(df)

  if (verbose) message("Reading annotation file: ", basename(annot_path))
  annot <- read.delim(annot_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
  .validate_lfq_annot(annot)
  annot$Column    <- as.character(annot$Column)
  annot$Condition <- as.character(annot$Condition)
  annot$Experiment <- if ("Experiment" %in% names(annot)) {
    as.character(annot$Experiment)
  } else {
    NA_character_
  }

  # --- Resolve the condition order and filter the design ---
  if (is.null(condition_order)) {
    condition_order <- unique(annot$Condition)
  } else {
    keep_ann <- annot$Condition %in% condition_order
    if (!any(keep_ann)) {
      stop("No condition in the _Annot matches condition_order = c(",
           paste0("'", condition_order, "'", collapse = ", "), ").\n",
           "Conditions in the _Annot: ", paste(unique(annot$Condition), collapse = ", "))
    }
    annot <- annot[keep_ann, , drop = FALSE]
  }

  # --- Sample table (sorted by condition and by the order in the Annot) ---
  annot$R.Condition <- factor(annot$Condition, levels = condition_order,
                              ordered = TRUE)
  annot <- annot[order(annot$R.Condition, seq_len(nrow(annot))), , drop = FALSE]
  # Replicate = sequential index within each condition (robust to names)
  annot$R.Replicate <- as.integer(
    stats::ave(seq_len(nrow(annot)), annot$R.Condition,
               FUN = function(i) seq_along(i))
  )
  coding_levels <- annot$Column   # Coding = sample name from the _Annot

  # --- Check that every sample in the design has its Abundance column ---
  abund_map <- .lfq_sample_cols(df, "^Abundance:\\s*", coding_levels)
  missing_abund <- setdiff(coding_levels, names(abund_map))
  if (length(missing_abund) > 0) {
    stop("No 'Abundance:' column found for these samples of the _Annot:\n  - ",
         paste(missing_abund, collapse = "\n  - "))
  }
  # Warn about Abundance columns in the TSV that are not part of the design
  all_abund <- sub("^Abundance:\\s*", "",
                   grep("^Abundance:\\s*", names(df), value = TRUE))
  extra_abund <- setdiff(all_abund, coding_levels)
  if (length(extra_abund) > 0 && verbose) {
    message("Note: 'Abundance:' columns ignored (not present in the _Annot): ",
            paste(extra_abund, collapse = ", "))
  }

  # ==========================================================================
  # metadata
  # ==========================================================================
  if (verbose) message("Generating sample metadata...")
  run_summary <- data.frame(
    R.FileName   = unname(abund_map[coding_levels]),
    R.Condition  = annot$R.Condition,
    R.Replicate  = annot$R.Replicate,
    Coding       = coding_levels,
    R.Experiment = annot$Experiment,
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
  col_qvalue <- .pd_resolve_col(df, "Exp. q-value: Combined", "Exp. q-value")
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

  abund_mat <- df[, unname(abund_map[coding_levels]), drop = FALSE]
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
    stop("Integrity error: protein_ID contains duplicated Accession values.")
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
    stop("Integrity error: protein_QUANT contains duplicated Accession values.")
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
  return(result)
}

# =============================================================================
# Methods for the lfq_data class
# =============================================================================

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
