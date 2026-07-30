# =============================================================================
# Spectronaut Data Preprocessing
# =============================================================================
#
# Converts Spectronaut reports into a structured format for downstream
# proteomics analysis.
#
# Copyright 2025 Sergio Ciordia
# License: GPL-3
# =============================================================================

# --- Dependencies ---

# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Get an aggregating function by name
#' @noRd
.get_aggregator <- function(name) {
  nm <- tolower(name %||% "")
  switch(nm,
    "max"    = function(x) max(x, na.rm = TRUE),
    "mean"   = function(x) mean(x, na.rm = TRUE),
    "median" = function(x) stats::median(x, na.rm = TRUE),
    "min"    = function(x) min(x, na.rm = TRUE),
    stop("Unsupported aggregator: '", name, "'. Use: max, mean, median, min.")
  )
}

#' Return the first non-NA value of a vector
#' @noRd
.first_non_na <- function(x) {

  y <- na.omit(x)
  if (length(y) == 0) NA else y[1]
}

#' Clean fields holding semicolon-separated values
#' @description Converts fields such as "16,8%;16.8%" to numeric and aggregates
#'   them by group
#' @noRd
.clean_semicolon_numeric <- function(df, value_col, group_cols, out_col = value_col,
                                     agg_fun = .get_aggregator("max")) {
  df %>%
    select(all_of(c(group_cols, value_col))) %>%
    mutate(!!value_col := as.character(.data[[value_col]])) %>%
    separate_rows(all_of(value_col), sep = ";") %>%
    mutate(
      .val_num = suppressWarnings(
        as.numeric(gsub("[^0-9.]+", "", gsub(",", ".", .data[[value_col]])))
      )
    ) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      !!out_col := if (all(is.na(.val_num))) NA_real_ else agg_fun(.val_num),
      .groups = "drop"
    )
}

#' Normalize Spectronaut column names
#' @description Standardises column-name variants to a consistent format
#' @noRd
.normalize_column_names <- function(df) {
  nm <- names(df)


  # PG.Cscore.RunWise
  cscore_variants <- c("PG.Cscore..Run.Wise.", "PG.Cscore (Run-Wise)", "PG.Cscore_RunWise")
  for (variant in cscore_variants) {
    if (variant %in% nm) {
      names(df)[nm == variant] <- "PG.Cscore.RunWise"
      nm <- names(df)
      break
    }
  }


  # Global / experiment-wide columns
  global_mappings <- list(
    "PG.NrOfPrecursorsIdentified..Experiment.wide." = "PG.NrOfPrecursorsIdentified.Global",
    "PG.NrOfStrippedSequencesIdentified..Experiment.wide." = "PG.NrOfStrippedSequencesIdentified.Global",
    "PG.Coverage..Global." = "PG.Coverage.Global"
  )

  for (old_name in names(global_mappings)) {
    if (old_name %in% nm) {
      names(df)[nm == old_name] <- global_mappings[[old_name]]
      nm <- names(df)
    }
  }

  df
}

#' Build the Coding column and sort by condition/replicate
#' @noRd
.make_coding <- function(df, cond_order) {
  df %>%
    mutate(
      R.Condition = factor(R.Condition, levels = cond_order, ordered = TRUE),
      R.Replicate = as.integer(R.Replicate),
      Coding = paste(R.Condition, R.Replicate, sep = "_")
    ) %>%
    arrange(R.Condition, R.Replicate)
}

#' Validate the required columns in a data frame
#' @noRd
.validate_spectronaut_columns <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(
      "Required columns missing from the file:\n",
      "  - ", paste(missing, collapse = "\n  - "), "\n",
      "Check that the file is a valid Spectronaut report."
    )
  }
  invisible(TRUE)
}

#' Extract the base protein information (static, 1 row per protein)
#' @noRd
.extract_base_info <- function(df, mw_clean) {
  static_cols <- c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight")

  df %>%
    select(all_of(static_cols)) %>%
    group_by(PG.ProteinGroups) %>%
    summarise(
      PG.ProteinDescriptions = .first_non_na(PG.ProteinDescriptions),
      PG.Genes = .first_non_na(PG.Genes),
      PG.MolecularWeight = .first_non_na(PG.MolecularWeight),
      .groups = "drop"
    ) %>%
    select(-PG.MolecularWeight) %>%
    left_join(mw_clean, by = "PG.ProteinGroups")
}

# =============================================================================
# Main Function
# =============================================================================

#' Preprocess Spectronaut reports
#'
#' @description
#' Converts a Spectronaut report (long TSV format) into three structured tables:
#' run metadata, protein identification and protein quantification.
#'
#' @param file_path Path to the Spectronaut TSV file.
#' @param condition_order Character vector with the order of the experimental
#'   conditions (e.g. `c("Control", "Treated")`).
#' @param export_dir Directory to export the TSV files to. If `NULL` (default),
#'   no files are exported.
#' @param agg_coverage_run Aggregation method for the per-sample PG.Coverage.
#'   Options: "max" (default), "mean", "median", "min".
#' @param agg_coverage_global Aggregation method for PG.Coverage.Global.
#'   Options: "max" (default), "mean", "median", "min".
#' @param agg_mw Aggregation method for PG.MolecularWeight.
#'   Options: "max" (default), "mean", "median", "min".
#' @param agg_cscore_runwise Aggregation method for PG.Cscore.RunWise.
#'   Options: "mean" (default), "max", "median", "min".
#' @param timestamp_suffix Logical. If `TRUE` (default), appends a timestamp to
#'   the names of the exported files.
#' @param verbose Logical. If `TRUE` (default), shows progress messages.
#'
#' @return A list with class `spectronaut_data` containing:
#'   \describe{
#'     \item{metadata}{Data frame with the run information (1 row per sample)}
#'     \item{protein_id}{Data frame with the identification metrics per protein}
#'     \item{protein_quant}{Data frame with the quantification metrics per
#'       protein}
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- preprocess_spectronaut(
#'   file_path = "data-raw/Spectronaut_Report.tsv",
#'   condition_order = c("Control", "Treatment")
#' )
#'
#' # Access the components
#' head(result$metadata)
#' head(result$protein_quant)
#'
#' # With export
#' result <- preprocess_spectronaut(
#'   file_path = "data-raw/Spectronaut_Report.tsv",
#'   condition_order = c("A", "B", "C", "D"),
#'   export_dir = "./results",
#'   agg_coverage_run = "mean",
#'   verbose = TRUE
#' )
#' }
#'
#' @export
preprocess_spectronaut <- function(
    file_path,
    condition_order,
    export_dir = NULL,
    agg_coverage_run = c("max", "mean", "median", "min"),
    agg_coverage_global = c("max", "mean", "median", "min"),
    agg_mw = c("max", "mean", "median", "min"),
    agg_cscore_runwise = c("mean", "max", "median", "min"),
    timestamp_suffix = TRUE,
    verbose = TRUE
) {

  # --- Argument validation ---
  agg_coverage_run <- match.arg(agg_coverage_run)
  agg_coverage_global <- match.arg(agg_coverage_global)
  agg_mw <- match.arg(agg_mw)
  agg_cscore_runwise <- match.arg(agg_cscore_runwise)

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  if (length(condition_order) == 0 || !is.character(condition_order)) {
    stop("condition_order must be a non-empty character vector.")
  }

  # --- Reading the file ---
  if (verbose) message("Reading file: ", basename(file_path))

  df <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = TRUE
  )

  # Normalize the column names

  df <- .normalize_column_names(df)

  # Validate the required columns
  required_cols <- c(
    "R.FileName", "R.Condition", "R.Replicate",
    "PG.ProteinGroups", "PG.Quantity"
  )
  .validate_spectronaut_columns(df, required_cols)

  # Create the Coding column
  df <- .make_coding(df, condition_order)

  # Get the sorted Coding levels (computed only once)
  coding_levels <- df %>%
    distinct(R.Condition, R.Replicate, Coding) %>%
    arrange(R.Condition, R.Replicate) %>%
    pull(Coding)

  # --- Create run_summary (metadata) ---
  if (verbose) message("Generating run metadata...")

  summary_cols <- c(
    "R.FileName", "R.Condition", "R.Replicate", "Coding",
    "R.PrecursorsIdentified", "R.StrippedSequencesIdentified",
    "R.ProteinGroupsIdentified"
  )

  run_summary <- df %>%
    select(any_of(summary_cols)) %>%
    distinct() %>%
    arrange(R.FileName) %>%
    group_by(R.FileName) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    arrange(R.Condition, R.Replicate) %>%
    as.data.frame()

  rownames(run_summary) <- run_summary$Coding

  # --- Prepare the aggregators ---
  cov_fun_run <- .get_aggregator(agg_coverage_run)
  cov_fun_global <- .get_aggregator(agg_coverage_global)
  mw_fun <- .get_aggregator(agg_mw)
  cscore_agg <- .get_aggregator(agg_cscore_runwise)

  # --- Cleaned MW (shared between protein_ID and protein_QUANT) ---
  mw_clean <- .clean_semicolon_numeric(
    df, "PG.MolecularWeight",
    group_cols = "PG.ProteinGroups",
    out_col = "PG.MolecularWeight",
    agg_fun = mw_fun
  )

  # ==========================================================================
  # protein_ID
  # ==========================================================================
  if (verbose) message("Processing protein_ID...")

  # Cleaned run-wise coverage

  coverage_clean <- .clean_semicolon_numeric(
    df, "PG.Coverage",
    group_cols = c("PG.ProteinGroups", "Coding"),
    out_col = "PG.Coverage",
    agg_fun = cov_fun_run
  )

  # Merge in the cleaned coverage and cast the numeric columns
  df2 <- df %>%
    select(-PG.Coverage) %>%
    left_join(coverage_clean, by = c("PG.ProteinGroups", "Coding")) %>%
    mutate(
      PG.NrOfPrecursorsIdentified = suppressWarnings(as.numeric(PG.NrOfPrecursorsIdentified)),
      PG.NrOfStrippedSequencesIdentified = suppressWarnings(as.numeric(PG.NrOfStrippedSequencesIdentified)),
      PG.Cscore.RunWise = suppressWarnings(as.numeric(PG.Cscore.RunWise))
    )

  # Base information
  base_info <- .extract_base_info(df2, mw_clean)

  # Run-wise metrics to wide format
  runwise_cols <- c(
    "PG.NrOfPrecursorsIdentified", "PG.NrOfStrippedSequencesIdentified",
    "PG.Coverage", "PG.Cscore.RunWise"
  )

  runwise_wide <- df2 %>%
    select(PG.ProteinGroups, Coding, all_of(runwise_cols)) %>%
    pivot_longer(cols = all_of(runwise_cols), names_to = "metric", values_to = "value") %>%
    group_by(PG.ProteinGroups, Coding, metric) %>%
    summarise(
      value = if (all(is.na(value))) NA_real_ else
        if (unique(metric) == "PG.Cscore.RunWise") cscore_agg(value) else
          if (unique(metric) == "PG.Coverage") cov_fun_run(value) else
            max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Coding = factor(Coding, levels = coding_levels, ordered = TRUE),
      metric_coding = paste0(metric, "_", Coding)
    ) %>%
    select(PG.ProteinGroups, metric_coding, value) %>%
    pivot_wider(names_from = metric_coding, values_from = value, names_repair = "unique")

  # Assemble protein_ID
  protein_ID <- base_info %>%
    left_join(runwise_wide, by = "PG.ProteinGroups")

  # Order the columns
  static_cols <- c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight")
  metric_order <- c(
    "PG.NrOfPrecursorsIdentified", "PG.NrOfStrippedSequencesIdentified",
    "PG.Coverage", "PG.Cscore.RunWise"
  )
  runwise_order <- unlist(lapply(metric_order, function(m) paste0(m, "_", coding_levels)))
  final_cols <- c(static_cols, intersect(runwise_order, names(protein_ID)))

  protein_ID <- protein_ID %>%
    select(all_of(final_cols)) %>%
    arrange(PG.ProteinGroups)

  # Check uniqueness
  if (n_distinct(protein_ID$PG.ProteinGroups) != nrow(protein_ID)) {
    stop(
      "Integrity error: protein_ID contains duplicated rows per PG.ProteinGroups. ",
      "Check the input data."
    )
  }

  # ==========================================================================
  # protein_QUANT
  # ==========================================================================
  if (verbose) message("Processing protein_QUANT...")

  # Cleaned Coverage.Global
  coverage_global_clean <- .clean_semicolon_numeric(
    df, "PG.Coverage.Global",
    group_cols = "PG.ProteinGroups",
    out_col = "PG.Coverage.Global",
    agg_fun = cov_fun_global
  )

  # Prepare the df for QUANT
  df4 <- df %>%
    select(-any_of("PG.Coverage.Global")) %>%
    left_join(coverage_global_clean, by = "PG.ProteinGroups") %>%
    mutate(
      PG.NrOfPrecursorsIdentified.Global = suppressWarnings(as.numeric(PG.NrOfPrecursorsIdentified.Global)),
      PG.NrOfStrippedSequencesIdentified.Global = suppressWarnings(as.numeric(PG.NrOfStrippedSequencesIdentified.Global)),
      PG.NrOfPrecursorsUsedForQuantification = suppressWarnings(as.numeric(PG.NrOfPrecursorsUsedForQuantification)),
      PG.NrOfStrippedSequencesUsedForQuantification = suppressWarnings(as.numeric(PG.NrOfStrippedSequencesUsedForQuantification)),
      PG.Quantity = suppressWarnings(as.numeric(PG.Quantity)),
      PG.Cscore = suppressWarnings(as.numeric(PG.Cscore))
    )

  # Base information (reuses mw_clean)
  base_info2 <- .extract_base_info(df4, mw_clean)

  # Global metrics
  global_metrics <- df4 %>%
    group_by(PG.ProteinGroups) %>%
    summarise(
      PG.NrOfPrecursorsIdentified.Global = .first_non_na(PG.NrOfPrecursorsIdentified.Global),
      PG.NrOfStrippedSequencesIdentified.Global = .first_non_na(PG.NrOfStrippedSequencesIdentified.Global),
      PG.Coverage.Global = .first_non_na(PG.Coverage.Global),
      PG.Cscore = .first_non_na(PG.Cscore),
      .groups = "drop"
    )

  # Per-sample metrics to wide format
  pivot_metrics <- c(
    "PG.NrOfPrecursorsUsedForQuantification",
    "PG.NrOfStrippedSequencesUsedForQuantification",
    "PG.Quantity"
  )

  runwise_wide2 <- df4 %>%
    select(PG.ProteinGroups, Coding, all_of(pivot_metrics)) %>%
    distinct() %>%
    pivot_longer(cols = all_of(pivot_metrics), names_to = "metric", values_to = "value") %>%
    group_by(PG.ProteinGroups, Coding, metric) %>%
    summarise(
      value = if (all(is.na(value))) NA_real_ else .first_non_na(value),
      .groups = "drop"
    ) %>%
    mutate(
      Coding = factor(Coding, levels = coding_levels, ordered = TRUE),
      metric_coding = paste0(metric, "_", Coding)
    ) %>%
    select(PG.ProteinGroups, metric_coding, value) %>%
    pivot_wider(names_from = metric_coding, values_from = value, names_repair = "unique")

  # Assemble protein_QUANT
  protein_QUANT <- base_info2 %>%
    left_join(global_metrics, by = "PG.ProteinGroups") %>%
    left_join(runwise_wide2, by = "PG.ProteinGroups")

  # Order the columns
  static_cols2 <- c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight")
  global_cols <- c(
    "PG.NrOfPrecursorsIdentified.Global",
    "PG.NrOfStrippedSequencesIdentified.Global",
    "PG.Coverage.Global",
    "PG.Cscore"
  )
  runwise_order2 <- unlist(lapply(pivot_metrics, function(m) paste0(m, "_", coding_levels)))
  final_cols2 <- c(static_cols2, global_cols, intersect(runwise_order2, names(protein_QUANT)))

  protein_QUANT <- protein_QUANT %>%
    select(all_of(final_cols2)) %>%
    arrange(PG.ProteinGroups)

  # Check uniqueness
  if (n_distinct(protein_QUANT$PG.ProteinGroups) != nrow(protein_QUANT)) {
    stop(
      "Integrity error: protein_QUANT contains duplicated rows per PG.ProteinGroups. ",
      "Check the input data."
    )
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
      "  - Runs: ", nrow(run_summary), "\n",
      "  - Proteins (ID): ", nrow(protein_ID), "\n",
      "  - Proteins (QUANT): ", nrow(protein_QUANT)
    )
  }

  result <- list(
    metadata = as.data.frame(run_summary),
    protein_id = as.data.frame(protein_ID),
    protein_quant = as.data.frame(protein_QUANT)
  )

  class(result) <- c("spectronaut_data", "proteomics_data", "list")
  return(result)
}

# =============================================================================
# Methods for the spectronaut_data class
# =============================================================================

#' @export
print.spectronaut_data <- function(x, ...) {
  cat("Preprocessed Spectronaut data\n")
  cat("-----------------------------\n")
  cat("Runs (metadata):", nrow(x$metadata), "\n")
  cat("Proteins (ID):", nrow(x$protein_id), "\n")
  cat("Proteins (QUANT):", nrow(x$protein_quant), "\n")
  cat("\nConditions:", paste(unique(x$metadata$R.Condition), collapse = ", "), "\n")
  invisible(x)
}
