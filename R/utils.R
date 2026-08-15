# =============================================================================
# Shared internal utilities
# =============================================================================
#
# Definitions common to every NADIA module. Each module used to redefine `%||%`
# on its own (18 copies with 3 different semantics), so the last definition
# loaded won and silently changed the behaviour of the modules already loaded.
# There is now a single canonical definition.
#
# The modules that need the RNG helpers loaded this file through the
# `if (!exists(".rng_state", mode = "function"))` block in their header. That
# block disappears once the project becomes a package: all the files under R/
# share one and the same namespace.
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================


#' Null-coalescing operator
#'
#' Returns `a` unless it is `NULL`, in which case it returns `b`. Same semantics
#' as `base::\%||\%` (available since R 4.4) and as `rlang::\%||\%`.
#'
#' It is defined here rather than imported from `base` so that the package keeps
#' working under the declared `Depends: R (>= 4.4)` without relying on the
#' operator being exported in that particular version. The definition is
#' identical.
#'
#' @param a Value to check.
#' @param b Fallback value if `a` is `NULL`.
#' @return `a` if it is not `NULL`; otherwise `b`.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Capture the state of the random number generator
#'
#' Functions that call `set.seed()` alter the RNG of the user's session, so any
#' code run afterwards is no longer reproducible. These two helpers make it
#' possible to leave the RNG exactly as it was:
#'
#' ```
#' old_rng <- .rng_state()
#' on.exit(.rng_restore(old_rng), add = TRUE)
#' set.seed(seed)
#' ```
#'
#' @return The value of `.Random.seed`, or `NULL` if the RNG has not been used
#'   yet in this session.
#' @keywords internal
#' @noRd
.rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

#' Restore the state of the random number generator
#'
#' @param state Value previously returned by `.rng_state()`.
#' @return `NULL`, invisibly. Called for its side effect.
#' @keywords internal
#' @noRd
.rng_restore <- function(state) {
  if (is.null(state)) {
    # The RNG had not been initialised before the call: the seed created by
    # set.seed() is removed so that no trace is left behind.
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}


# =============================================================================
# Colour helpers
# =============================================================================
#
# These were duplicated in Boxplot_Highcharts_Final.R, PCA_Highcharts_Final.R,
# Pattern_Profiler_Highcharts.R and Heatmap_tidyHeatmap.R (8 definitions of 3
# functions). With every module living in the global environment, the active
# copy was the one from the last module loaded; in a package the winner would
# be the one from the last file in alphabetical order. They are unified here.
#
# Verified before unifying: for hexadecimal inputs the three copies of
# `hex_to_rgba` and the two of `darken_hex` produced identical output for the
# same alpha/factor. They only differed in the default value, and all 7 call
# sites always pass it explicitly. That is why NO default is declared here: a
# call without the argument must fail visibly instead of taking some arbitrary
# value.


#' Normalise a hexadecimal colour
#'
#' Strips the leading `#`, discards the alpha channel if the colour comes in
#' `RRGGBBAA` format, and returns `#RRGGBB` in upper case.
#'
#' @param hex Colour in hexadecimal format, with or without `#`.
#' @return `#RRGGBB` string in upper case.
#' @keywords internal
#' @noRd
.normalize_hex <- function(hex) {
  hex <- gsub("^#", "", hex)
  if (nchar(hex) == 8) {
    hex <- substr(hex, 1, 6)
  }
  paste0("#", toupper(hex))
}


#' Convert a hexadecimal colour to an `rgba()` string
#'
#' @param hex Hexadecimal colour.
#' @param alpha Opacity between 0 and 1. Deliberately without a default value:
#'   every call site specifies it.
#' @return `"rgba(r, g, b, a)"` string suitable for Highcharts.
#' @keywords internal
#' @noRd
.hex_to_rgba <- function(hex, alpha) {
  hex <- .normalize_hex(hex)
  rgb_vals <- grDevices::col2rgb(hex)
  sprintf("rgba(%d, %d, %d, %.2f)",
          rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha)
}


#' Darken a hexadecimal colour
#'
#' @param hex Hexadecimal colour.
#' @param factor Darkening fraction between 0 and 1. Deliberately without a
#'   default value: every call site specifies it.
#' @return Darkened `#RRGGBB` colour.
#' @keywords internal
#' @noRd
.darken_hex <- function(hex, factor) {
  hex <- .normalize_hex(hex)
  rgb_vals <- grDevices::col2rgb(hex)
  rgb_dark <- pmax(0, rgb_vals * (1 - factor))
  sprintf("#%02X%02X%02X",
          round(rgb_dark[1]), round(rgb_dark[2]), round(rgb_dark[3]))
}


# =============================================================================
# Feature filtering helpers
# =============================================================================

#' Build the adjusted p-value column name for a comparison
#'
#' @param comparison Name of the comparison, e.g. `"B-A"`.
#' @return Column name, e.g. `"adjP_B-A"`.
#' @keywords internal
#' @noRd
.adjp_col <- function(comparison) {
  paste0("adjP_", comparison)
}


#' Get the feature IDs according to the filtering mode
#'
#' Unifies the two copies that lived in `Heatmap_tidyHeatmap.R` and
#' `PCA_Highcharts_Final.R`. They differed in two respects:
#'
#' * The name of the third mode: `"target"` in the heatmap copy and `"specific"`
#'   in the PCA one. **Both** are accepted here as synonyms, because the public
#'   functions of each module propagate their own vocabulary and either one
#'   would have broken the other when merged into a single namespace.
#' * The `mode = "any"` filtering: the PCA copy used `sig_any == TRUE`, which
#'   lets `NA` `FeatureID`s slip through when `sig_any` contains `NA`. The
#'   `which()` from the heatmap copy is kept, which does not.
#'
#' @param data Data frame in long format with columns `FeatureID`, `sig_any` and
#'   `adjP_*`.
#' @param mode `"all"` (everything), `"any"` (significant in at least one
#'   comparison) or `"target"`/`"specific"` (significant in `comparison`).
#' @param alpha Significance threshold for the targeted mode.
#' @param comparison Comparison to use in the targeted mode.
#' @return Vector of `FeatureID`s meeting the criterion.
#' @keywords internal
#' @noRd
.get_feature_ids <- function(data,
                             mode = c("all", "any", "target", "specific"),
                             alpha = 0.05,
                             comparison = NULL) {

  mode <- match.arg(mode)
  data <- as.data.frame(data)

  feat <- data[!duplicated(data$FeatureID), , drop = FALSE]

  if (mode == "all") {
    return(feat$FeatureID)
  }

  if (mode == "any") {
    if (!("sig_any" %in% names(feat))) {
      stop("Column 'sig_any' is required for mode = 'any'")
    }
    # which() prevents NA FeatureIDs from slipping through when sig_any contains
    # NA (unlike feat$sig_any == TRUE, which would return NA rows).
    return(feat$FeatureID[which(feat$sig_any)])
  }

  # targeted mode: "target" and "specific" are synonyms
  if (is.null(comparison)) {
    stop("Argument 'comparison' is required for mode = '", mode, "'")
  }

  col <- .adjp_col(comparison)
  if (!(col %in% names(feat))) {
    stop("Column not found: ", col)
  }

  feat$FeatureID[which(feat[[col]] <= alpha)]
}


# =============================================================================
# Mfuzz dependencies that must be attached
# =============================================================================
#
# `Mfuzz` declares `Depends: Biobase, e1071` and calls `exprs()` and `cmeans()`
# unqualified, so its functions only resolve those names if both packages are on
# the search path. `requireNamespace()` is not enough: it loads the namespace but
# does not attach it. This used to be taken care of by the `library(Mfuzz)` call
# in the module header, which pulled in its Depends; a package cannot do that.
#
# The solution is to attach them only while the Pattern Profiler is running and
# to leave the search path as it was, so that the user's session is untouched.

#' Attach the dependencies that Mfuzz needs on the search path
#'
#' @return Vector with the packages that this call has attached (possibly empty
#'   if they were already attached). It must be passed to
#'   `.mfuzz_deps_detach()`.
#' @keywords internal
#' @noRd
.mfuzz_deps_attach <- function() {
  required <- c("Mfuzz", "Biobase", "e1071")
  missing_pkgs <- required[!vapply(required, requireNamespace, logical(1),
                                   quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("The Pattern Profiler requires ", paste(missing_pkgs, collapse = ", "),
         ". Install them with BiocManager::install(c(",
         paste(sprintf('"%s"', missing_pkgs), collapse = ", "), ")).",
         call. = FALSE)
  }
  # Only Biobase and e1071 need to be attached; Mfuzz is used via Mfuzz::
  attachable <- c("Biobase", "e1071")
  already_on <- paste0("package:", attachable) %in% search()
  for (p in attachable[!already_on]) attachNamespace(asNamespace(p))
  attachable[!already_on]
}

#' Undo what `.mfuzz_deps_attach()` did
#'
#' @param pkgs Vector returned by `.mfuzz_deps_attach()`.
#' @return `NULL`, invisibly. Called for its side effect.
#' @keywords internal
#' @noRd
.mfuzz_deps_detach <- function(pkgs) {
  for (p in pkgs) {
    entry <- paste0("package:", p)
    if (entry %in% search()) {
      detach(entry, character.only = TRUE, unload = FALSE)
    }
  }
  invisible(NULL)
}


# =============================================================================
# Shared preprocessing helpers
# =============================================================================

#' Detect the `Abundance:` columns and derive Coding / Condition / Replicate
#'
#' Two of the wide formats declare their design in the column names themselves:
#' Proteome Discoverer writes `Abundance: <Condition>_<Replicate>`, and a DIA-NN
#' matrix can be renamed to the same convention so that it needs no annotation
#' file either. The parser is therefore shared by `preprocess_tmt()` and
#' `preprocess_diann()`.
#'
#' @param df Data frame read from the export.
#' @return A data frame with `abundance_col`, `Coding`, `R.Condition` and
#'   `R.Replicate`, one row per column found.
#' @keywords internal
#' @noRd
.parse_abundance_columns <- function(df) {
  abund_cols <- grep("^Abundance:\\s*", names(df), value = TRUE)
  if (length(abund_cols) == 0) {
    stop("No 'Abundance:' columns found in the file. ",
         "Check that it is a Proteome Discoverer export.")
  }

  coding <- sub("^Abundance:\\s*", "", abund_cols)
  m <- stringr::str_match(coding, "^(.+)_(\\d+)$")
  if (any(is.na(m[, 1]))) {
    bad <- coding[is.na(m[, 1])]
    stop(
      "Could not parse the Abundance suffixes (expected <Condition>_<Replicate>):\n  - ",
      paste(bad, collapse = "\n  - ")
    )
  }

  data.frame(
    abundance_col = abund_cols,
    Coding        = coding,
    R.Condition   = m[, 2],
    R.Replicate   = as.integer(m[, 3]),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Shared by the four preprocessing functions
# =============================================================================

#' Drop conditions that were asked for but are not in the data
#'
#' @description
#' Every `preprocess_*()` function takes `condition_order`, uses it to select the
#' samples to keep, and then makes it the levels of the `R.Condition` factor.
#' Those two roles disagree when a condition is listed but does not occur in the
#' file: the selection finds nothing to keep, while the factor gains a level with
#' no samples behind it. An empty level is not inert -- it travels into
#' `colData()` and from there into legends and axes.
#'
#' Rather than fail (the other conditions are perfectly usable), this warns and
#' narrows `condition_order` to what is actually there. `setdiff()` preserves the
#' order of its first argument, so the order that was asked for survives.
#'
#' @param condition_order Character vector as supplied by the caller.
#' @param present The conditions occurring in the data, in any form.
#' @param source Word naming where they were looked for, for the message.
#' @return `condition_order` without the absent conditions.
#' @keywords internal
#' @noRd
.drop_absent_conditions <- function(condition_order, present, source = "report") {
  absent <- setdiff(condition_order, unique(as.character(present)))
  if (length(absent) > 0) {
    warning("Conditions listed in condition_order but absent from the ", source,
            ", and dropped: ", paste(absent, collapse = ", "), call. = FALSE)
    condition_order <- setdiff(condition_order, absent)
  }
  condition_order
}

#' Read a sample sheet and turn it into a design table
#'
#' @description
#' The escape hatch for reports whose intensity columns are not named
#' `<Condition>_<Replicate>`. It is a file rather than a vector of labels on
#' purpose: a vector is positional, and column order is not something an export
#' guarantees, so a re-export that reorders the columns would relabel every run
#' without any way to notice. Matching by name either matches or fails.
#'
#' Replicate numbers come from the identifier when it ends in `_<digits>`, which
#' is what the suffix route does, so a sheet stating the obvious gives exactly
#' the same answer as no sheet at all. Identifiers that do not carry a number --
#' raw-file names, TMT tags -- are numbered sequentially within their condition.
#' Either way the number does not depend on the order the sheet lists them in,
#' so a reordered sheet cannot relabel anything.
#'
#' @param annot_path Path to a TSV with `Column` and `Condition`. Any other
#'   column is ignored: batch and covariate structure belongs in `covariate_df`
#'   of [process_proteomics()], not here.
#' @param ids Sample identifiers as they appear in the report, one per intensity
#'   column, in column order.
#' @param cols The intensity column names themselves, parallel to `ids`.
#' @param verbose Report the columns the sheet leaves out.
#' @return A data frame with `abundance_col`, `Coding`, `R.Condition` and
#'   `R.Replicate` -- the same shape `.parse_abundance_columns()` returns, so
#'   that everything downstream is blind to which route was taken.
#' @keywords internal
#' @noRd
.design_from_annotation <- function(annot_path, ids, cols, verbose = TRUE) {
  if (!file.exists(annot_path)) stop("Annotation file not found: ", annot_path)

  annot <- utils::read.delim(annot_path, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE)

  missing_cols <- setdiff(c("Column", "Condition"), names(annot))
  if (length(missing_cols) > 0) {
    stop("Required columns missing from the annotation file:\n  - ",
         paste(missing_cols, collapse = "\n  - "),
         "\nIt must have at least 'Column' and 'Condition'.")
  }
  annot$Column    <- as.character(annot$Column)
  annot$Condition <- as.character(annot$Condition)

  if (anyDuplicated(annot$Column) > 0) {
    stop("The annotation file has duplicated values in 'Column':\n  - ",
         paste(unique(annot$Column[duplicated(annot$Column)]),
               collapse = "\n  - "))
  }

  hit <- match(annot$Column, ids)
  if (anyNA(hit)) {
    stop("These samples of the annotation file have no column in the report:",
         "\n  - ", paste(annot$Column[is.na(hit)], collapse = "\n  - "),
         "\nThe report has:\n  - ", paste(ids, collapse = "\n  - "))
  }

  # Columns the sheet leaves out are dropped, which is how a pool or a blank is
  # excluded. Worth saying out loud, since it silently shrinks the experiment.
  extra <- setdiff(ids, annot$Column)
  if (length(extra) > 0 && verbose) {
    message("Note: columns ignored (not listed in the annotation file): ",
            paste(extra, collapse = ", "))
  }

  # Take the replicate from the identifier when it has one, so that a sheet
  # spelling out what the suffixes already say lands on the same answer. Fall
  # back to a per-condition counter for identifiers that carry no number, and
  # take that fallback for the whole sheet rather than per row, so the numbering
  # cannot be a mixture of two schemes.
  m <- stringr::str_match(annot$Column, "^(.+)_(\\d+)$")
  replicate <- if (anyNA(m[, 1])) {
    as.integer(stats::ave(seq_len(nrow(annot)), annot$Condition, FUN = seq_along))
  } else {
    as.integer(m[, 3])
  }

  data.frame(
    abundance_col = cols[hit],
    Coding        = annot$Column,
    R.Condition   = annot$Condition,
    R.Replicate   = replicate,
    stringsAsFactors = FALSE
  )
}

#' Resolve the design of a Proteome Discoverer export
#'
#' @description Suffixes if they are there, sample sheet if they are not. Shared
#'   by `preprocess_tmt()` and `preprocess_lfq()`, whose exports both name their
#'   intensity columns `Abundance: <sample>`.
#' @keywords internal
#' @noRd
.pd_resolve_design <- function(df, annot_path = NULL, verbose = TRUE) {
  if (is.null(annot_path)) return(.parse_abundance_columns(df))

  cols <- grep("^Abundance:\\s*", names(df), value = TRUE)
  if (length(cols) == 0) {
    stop("No 'Abundance:' columns found in the file. ",
         "Check that it is a Proteome Discoverer export.")
  }
  .design_from_annotation(annot_path, sub("^Abundance:\\s*", "", cols), cols,
                          verbose = verbose)
}

# =============================================================================
# Proteome Discoverer helpers
# =============================================================================
#
# These were duplicated between Preprocessing_LFQ.R and Preprocessing_TMT.R. The
# two copies of `.parse_gene_from_description` were identical; those of
# `.validate_pd_columns` differed only in the text of the error message.

#' Extract the Gene Name from the Description column (embedded UniProt format)
#'
#' Looks for the `GN=<gene>` pattern usual in Proteome Discoverer exports.
#'
#' @param x Vector of descriptions.
#' @return Vector of gene names, `NA` where there is no match.
#' @keywords internal
#' @noRd
.parse_gene_from_description <- function(x) {
  m <- stringr::str_match(x, "GN=([^ ]+)")
  m[, 2]
}

#' Record where an object came from
#'
#' Attaches the call that produced an object and a fingerprint of the files it
#' was read from, as the attributes `nadia_call` and `nadia_source`.
#'
#' They are attributes rather than list elements on purpose: every preprocessing
#' function documents a return value of exactly three elements (`metadata`,
#' `protein_id`, `protein_quant`), and adding a fourth would change a contract
#' that code and documentation already depend on. An attribute travels with the
#' object without appearing in `names()` or in `str()`.
#'
#' The fingerprint is what lets someone check, years later, whether a result
#' corresponds to the raw file in front of them. The MD5 is skipped for files
#' over 500 MB, where hashing costs more than the answer is worth.
#'
#' @param x     Object to stamp.
#' @param call  The call to record, normally `match.call()` from the caller.
#' @param files Named character vector of paths; the names become the `role`
#'   column (`report`, `annotation`).
#' @return `x`, with the two attributes attached.
#' @keywords internal
#' @noRd
.nadia_stamp <- function(x, call, files = character(0)) {
  attr(x, "nadia_call") <- call

  files <- files[!vapply(files, is.null, logical(1))]
  files <- unlist(files)

  if (length(files) > 0) {
    info <- file.info(files)
    attr(x, "nadia_source") <- data.frame(
      role       = names(files),
      path       = normalizePath(files, mustWork = FALSE),
      size_bytes = as.numeric(info$size),
      mtime      = format(info$mtime, tz = "UTC", usetz = FALSE),
      md5        = vapply(seq_along(files), function(i) {
        if (is.na(info$size[i]) || info$size[i] > 500e6) {
          NA_character_
        } else {
          unname(tools::md5sum(files[i]))
        }
      }, character(1)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  x
}


#' Validate the minimum columns of a Proteome Discoverer export
#'
#' @param df Data frame read from the export.
#' @return `TRUE` invisibly; aborts with an error if any required column is
#'   missing.
#' @keywords internal
#' @noRd
.validate_pd_columns <- function(df) {
  required <- c("Accession", "Description")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "Required columns missing from the data file:\n  - ",
      paste(missing, collapse = "\n  - "),
      "\nCheck that the file is a Proteome Discoverer protein export."
    )
  }
  invisible(TRUE)
}
