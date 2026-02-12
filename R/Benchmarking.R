# =============================================================================
# Benchmarking Module
# =============================================================================
#
# Benchmark differential expression results against known ground truth
# (multi-species spike-in experiments). Computes classification metrics
# (Sensitivity, Specificity, AUC, etc.), dispersion statistics, and
# generates interactive (Highcharter) and static (ggplot2) visualizations.
#
# Input:
#   - de_res: Data frame with columns Protein.IDs, Gene.Names, logFC,
#             P.Value, adj.P.Val, Change, Comparison, Assay, Species
#   - expected_values: Data frame with columns Comparison, Species,
#                      expected_logFC
#
# Dependencies: dplyr, tidyr, highcharter, ggplot2, scales
# Optional: pROC (for AUC), arrow (parquet), readr (TSV)
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# --- Null coalescing operator ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}


# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Validate required columns in DE results
#'
#' @param de_res Data frame with differential expression results
#' @return Invisible TRUE if valid, stops with error otherwise
#' @keywords internal
.validate_de_res <- function(de_res) {
  required <- c("Protein.IDs", "logFC", "Comparison", "Species")
  missing <- setdiff(required, names(de_res))
  if (length(missing) > 0) {
    stop("Columnas requeridas faltantes en de_res: ",
         paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

#' Validate required columns in expected_values
#'
#' @param ev Data frame with expected values
#' @return Invisible TRUE if valid, stops with error otherwise
#' @keywords internal
.validate_expected_values <- function(ev) {
  required <- c("Comparison", "Species", "expected_logFC")
  missing <- setdiff(required, names(ev))
  if (length(missing) > 0) {
    stop("Columnas requeridas faltantes en expected_values: ",
         paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

#' Auto-detect p-value column
#'
#' @param df Data frame
#' @param preferred Preferred column name
#' @return Name of the p-value column found
#' @keywords internal
.detect_p_col <- function(df, preferred = "adj.P.Val") {
  if (preferred %in% names(df)) return(preferred)
  candidates <- c("adj.P.Val", "P.Value", "pvalue", "p.value", "padj")
  found <- intersect(candidates, names(df))
  if (length(found) == 0) {
    stop("No se encontro columna de p-valor. ",
         "Columnas disponibles: ", paste(names(df), collapse = ", "))
  }
  found[1]
}

#' Harmonize Change column values
#'
#' Normalizes "Up Regulated" -> "Up", "Down Regulated" -> "Down",
#' "Not Significant" / "No Change" -> "No Change"
#'
#' @param de_res Data frame with Change column
#' @return Data frame with harmonized Change column
#' @keywords internal
.harmonize_change <- function(de_res) {
  if (!"Change" %in% names(de_res)) return(de_res)

  change <- as.character(de_res$Change)
  change[grepl("^Up", change, ignore.case = TRUE)] <- "Up"
  change[grepl("^Down", change, ignore.case = TRUE)] <- "Down"
  change[grepl("No Change|Not Significant|NS", change, ignore.case = TRUE)] <- "No Change"
  de_res$Change <- change
  de_res
}

#' Configure species colors
#'
#' @param species Character vector of unique species names
#' @param colors Named vector of custom colors (optional)
#' @return Named vector of colors per species
#' @keywords internal
.configure_species_colors <- function(species, colors = NULL) {
  defaults <- c(
    "ECOLI"  = "#D55E00",
    "YEAST"  = "#0072B2",
    "HUMAN"  = "#666666",
    "UPS"    = "#009E73",
    "ARATH"  = "#CC79A7",
    "BOVIN"  = "#E69F00"
  )

  result <- character(length(species))
  names(result) <- species

  for (sp in species) {
    if (!is.null(colors) && sp %in% names(colors)) {
      result[sp] <- colors[sp]
    } else if (sp %in% names(defaults)) {
      result[sp] <- defaults[sp]
    } else {
      fallback <- c("#56B4E9", "#F0E442", "#999999", "#882255",
                    "#332288", "#117733", "#44AA99", "#88CCEE")
      idx <- which(species == sp)
      result[sp] <- fallback[((idx - 1) %% length(fallback)) + 1]
    }
  }

  result
}

#' Convert hex color to rgba string
#'
#' @param hex Hexadecimal color string
#' @param alpha Opacity value (0-1)
#' @return String in rgba() format
#' @keywords internal
.hex_to_rgba <- function(hex, alpha = 1) {
  hex <- gsub("^#", "", hex)
  if (nchar(hex) == 8) hex <- substr(hex, 1, 6)
  hex <- paste0("#", toupper(hex))
  rgb_vals <- col2rgb(hex)
  sprintf("rgba(%d, %d, %d, %.2f)",
          rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha)
}

#' Compute trimmed SD and CV
#'
#' @param x Numeric vector
#' @param trim Proportion to trim from each side (default: 0.1)
#' @return Named list with trimmed_sd and trimmed_cv
#' @keywords internal
.trimmed_sd_cv <- function(x, trim = 0.1) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(list(trimmed_sd = NA_real_, trimmed_cv = NA_real_))

  n <- length(x)
  lo <- floor(n * trim) + 1
  hi <- n - floor(n * trim)
  x_sorted <- sort(x)
  x_trimmed <- x_sorted[lo:hi]

  tsd <- sd(x_trimmed)
  tmean <- mean(x_trimmed)
  tcv <- if (abs(tmean) > .Machine$double.eps) abs(tsd / tmean) * 100 else NA_real_

  list(trimmed_sd = tsd, trimmed_cv = tcv)
}


# =============================================================================
# DATA PREPARATION
# =============================================================================

#' Prepare benchmark data
#'
#' Filters de_res by assay and comparisons, validates inputs,
#' resolves p-value column.
#'
#' @param de_res Data frame with DE results (must include Species column)
#' @param expected_values Data frame with expected values
#' @param alpha Significance threshold
#' @param lfc_thr Log fold-change threshold
#' @param p_col P-value column name
#' @param comparisons Comparisons to include (NULL = all)
#' @param assay Assay to filter (NULL = all)
#' @return List with filtered de_res, expected_values, and resolved p_col
#' @keywords internal
.prepare_benchmark_data <- function(de_res, expected_values,
                                    alpha = 0.05, lfc_thr = 0,
                                    p_col = "adj.P.Val",
                                    comparisons = NULL, assay = NULL) {
  # Validate
  .validate_de_res(de_res)
  .validate_expected_values(expected_values)

  # Filter by assay
  if (!is.null(assay) && "Assay" %in% names(de_res)) {
    de_res <- de_res[de_res$Assay %in% assay, , drop = FALSE]
  }

  # Filter by comparisons
  if (!is.null(comparisons)) {
    de_res <- de_res[de_res$Comparison %in% comparisons, , drop = FALSE]
    expected_values <- expected_values[expected_values$Comparison %in% comparisons, , drop = FALSE]
  }

  if (nrow(de_res) == 0) stop("de_res vacio tras aplicar filtros.")
  if (nrow(expected_values) == 0) stop("expected_values vacio tras aplicar filtros.")

  # Harmonize Change column
  de_res <- .harmonize_change(de_res)

  # Resolve p-value column
  p_col <- .detect_p_col(de_res, p_col)

  # Ensure comparisons align
  de_comps <- unique(de_res$Comparison)
  ev_comps <- unique(expected_values$Comparison)
  common_comps <- intersect(de_comps, ev_comps)

  if (length(common_comps) == 0) {
    stop("No hay comparaciones en comun entre de_res y expected_values.\n",
         "  de_res: ", paste(de_comps, collapse = ", "), "\n",
         "  expected_values: ", paste(ev_comps, collapse = ", "))
  }

  de_res <- de_res[de_res$Comparison %in% common_comps, , drop = FALSE]
  expected_values <- expected_values[expected_values$Comparison %in% common_comps, , drop = FALSE]

  list(
    de_res = de_res,
    expected_values = expected_values,
    p_col = p_col,
    comparisons = common_comps
  )
}


# =============================================================================
# GROUND TRUTH CLASSIFICATION
# =============================================================================

#' Classify proteins for a single comparison
#'
#' Assigns truth labels and predicted labels based on expected values.
#' - Species in expected_values -> truth = 1 (expected change)
#' - Species NOT in expected_values -> truth = 0 (no expected change)
#' - predicted = 1 if significant AND direction matches expected_logFC sign
#' - For negative species: predicted = 1 if significant (any direction = FP)
#'
#' @param de_res_comp DE results for one comparison
#' @param ev_comp Expected values for one comparison
#' @param alpha Significance threshold
#' @param lfc_thr Log fold-change threshold
#' @param p_col P-value column name
#' @return Data frame with added columns: truth, predicted, classification
#' @keywords internal
.classify_proteins <- function(de_res_comp, ev_comp, alpha, lfc_thr, p_col) {
  # Species with expected changes (positives)
  positive_species <- ev_comp$Species
  expected_lfc <- setNames(ev_comp$expected_logFC, ev_comp$Species)

  # All species in data
  all_species <- unique(de_res_comp$Species)

  # Assign truth: 1 = expected change, 0 = no expected change
  de_res_comp$truth <- ifelse(de_res_comp$Species %in% positive_species, 1L, 0L)

  # Assign predicted
  pvals <- suppressWarnings(as.numeric(de_res_comp[[p_col]]))
  lfc <- as.numeric(de_res_comp$logFC)
  is_significant <- !is.na(pvals) & pvals <= alpha & abs(lfc) >= lfc_thr

  # For positive species: must have correct direction
  correct_direction <- logical(nrow(de_res_comp))
  for (sp in positive_species) {
    mask <- de_res_comp$Species == sp
    exp_sign <- sign(expected_lfc[sp])
    correct_direction[mask] <- sign(lfc[mask]) == exp_sign
  }
  # For negative species: any significant result counts as positive prediction
  negative_mask <- !(de_res_comp$Species %in% positive_species)

  de_res_comp$predicted <- 0L
  # Positive species: significant AND correct direction
  pos_mask <- de_res_comp$Species %in% positive_species
  de_res_comp$predicted[pos_mask & is_significant & correct_direction] <- 1L
  # Negative species: any significant = predicted positive (potential FP)
  de_res_comp$predicted[negative_mask & is_significant] <- 1L

  # Classification: TP, FP, TN, FN
  de_res_comp$classification <- "TN"
  de_res_comp$classification[de_res_comp$truth == 1 & de_res_comp$predicted == 1] <- "TP"
  de_res_comp$classification[de_res_comp$truth == 0 & de_res_comp$predicted == 1] <- "FP"
  de_res_comp$classification[de_res_comp$truth == 1 & de_res_comp$predicted == 0] <- "FN"
  de_res_comp$classification[de_res_comp$truth == 0 & de_res_comp$predicted == 0] <- "TN"

  de_res_comp
}

#' Classify proteins for all comparisons
#'
#' @param de_res Data frame with DE results
#' @param ev Data frame with expected values
#' @param alpha Significance threshold
#' @param lfc_thr Log fold-change threshold
#' @param p_col P-value column name
#' @return Data frame with classification columns for all comparisons
#' @keywords internal
.classify_all_comparisons <- function(de_res, ev, alpha, lfc_thr, p_col) {
  comparisons <- unique(ev$Comparison)

  classified_list <- lapply(comparisons, function(comp) {
    de_comp <- de_res[de_res$Comparison == comp, , drop = FALSE]
    ev_comp <- ev[ev$Comparison == comp, , drop = FALSE]

    if (nrow(de_comp) == 0 || nrow(ev_comp) == 0) return(NULL)

    .classify_proteins(de_comp, ev_comp, alpha, lfc_thr, p_col)
  })

  do.call(rbind, Filter(Negate(is.null), classified_list))
}


# =============================================================================
# METRICS COMPUTATION
# =============================================================================

#' Compute classification metrics from confusion matrix counts
#'
#' @param tp True positives
#' @param fp False positives
#' @param tn True negatives
#' @param fn False negatives
#' @return Named numeric vector with metrics
#' @keywords internal
.compute_metrics <- function(tp, fp, tn, fn) {
  # Convert to double to avoid integer overflow in MCC computation
  tp <- as.double(tp)
  fp <- as.double(fp)
  tn <- as.double(tn)
  fn <- as.double(fn)

  total <- tp + fp + tn + fn

  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  npv         <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  accuracy    <- if (total > 0) (tp + tn) / total else NA_real_
  f1          <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else {
    NA_real_
  }

  # Matthews Correlation Coefficient
  mcc_num <- (tp * tn) - (fp * fn)
  mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (mcc_den > 0) mcc_num / mcc_den else NA_real_

  c(
    TP          = tp,
    FP          = fp,
    TN          = tn,
    FN          = fn,
    Sensitivity = round(sensitivity, 4),
    Specificity = round(specificity, 4),
    Precision   = round(precision, 4),
    NPV         = round(npv, 4),
    Accuracy    = round(accuracy, 4),
    F1          = round(f1, 4),
    MCC         = round(mcc, 4)
  )
}

#' Compute AUC using pROC (optional)
#'
#' @param classified_df Classified data frame for one comparison
#' @param p_col P-value column name
#' @return AUC value or NA if pROC not available
#' @keywords internal
.compute_auc <- function(classified_df, p_col) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    return(NA_real_)
  }

  truth <- classified_df$truth
  predictor <- suppressWarnings(as.numeric(classified_df[[p_col]]))

  # Remove NAs
  valid <- !is.na(truth) & !is.na(predictor)
  truth <- truth[valid]
  predictor <- predictor[valid]

  if (length(unique(truth)) < 2) return(NA_real_)
  if (length(truth) < 3) return(NA_real_)

  tryCatch({
    roc_obj <- pROC::roc(truth, predictor, direction = ">", quiet = TRUE)
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) {
    NA_real_
  })
}

#' Compute benchmark metrics for all comparisons
#'
#' @param de_res Data frame with DE results (must include Species)
#' @param ev Data frame with expected values
#' @param alpha Significance threshold (default: 0.05)
#' @param lfc_thr Log fold-change threshold (default: 0)
#' @param p_col P-value column name (default: "adj.P.Val")
#' @param comparisons Comparisons to include (NULL = all)
#' @param assay Assay to filter (NULL = all)
#'
#' @return Data frame with metrics per comparison
#' @export
compute_benchmark_metrics <- function(de_res, ev,
                                      alpha = 0.05,
                                      lfc_thr = 0,
                                      p_col = "adj.P.Val",
                                      comparisons = NULL,
                                      assay = NULL) {
  # Prepare data
  prep <- .prepare_benchmark_data(de_res, ev, alpha, lfc_thr, p_col,
                                  comparisons, assay)
  de_res <- prep$de_res
  ev <- prep$expected_values
  p_col <- prep$p_col

  # Classify
  classified <- .classify_all_comparisons(de_res, ev, alpha, lfc_thr, p_col)

  comps <- unique(classified$Comparison)

  metrics_list <- lapply(comps, function(comp) {
    df_comp <- classified[classified$Comparison == comp, , drop = FALSE]

    tp <- sum(df_comp$classification == "TP")
    fp <- sum(df_comp$classification == "FP")
    tn <- sum(df_comp$classification == "TN")
    fn <- sum(df_comp$classification == "FN")

    m <- .compute_metrics(tp, fp, tn, fn)
    auc <- .compute_auc(df_comp, p_col)

    data.frame(
      Comparison  = comp,
      TP          = m["TP"],
      FP          = m["FP"],
      TN          = m["TN"],
      FN          = m["FN"],
      Sensitivity = m["Sensitivity"],
      Specificity = m["Specificity"],
      Precision   = m["Precision"],
      NPV         = m["NPV"],
      Accuracy    = m["Accuracy"],
      F1          = m["F1"],
      MCC         = m["MCC"],
      AUC         = round(auc, 4),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  do.call(rbind, metrics_list)
}

#' Compute confusion matrix by species
#'
#' @param classified_df Classified data frame
#' @return Data frame with TP/FP/TN/FN counts per Comparison x Species
#' @keywords internal
.confusion_by_species <- function(classified_df) {
  comps <- unique(classified_df$Comparison)
  species_all <- unique(classified_df$Species)

  result_list <- lapply(comps, function(comp) {
    df_comp <- classified_df[classified_df$Comparison == comp, , drop = FALSE]

    lapply(species_all, function(sp) {
      df_sp <- df_comp[df_comp$Species == sp, , drop = FALSE]
      if (nrow(df_sp) == 0) return(NULL)

      n_total <- nrow(df_sp)
      tp <- sum(df_sp$classification == "TP")
      fp <- sum(df_sp$classification == "FP")
      tn <- sum(df_sp$classification == "TN")
      fn <- sum(df_sp$classification == "FN")

      data.frame(
        Comparison = comp,
        Species    = sp,
        N          = n_total,
        TP         = tp,
        FP         = fp,
        TN         = tn,
        FN         = fn,
        TP_pct     = round(tp / n_total * 100, 1),
        FP_pct     = round(fp / n_total * 100, 1),
        TN_pct     = round(tn / n_total * 100, 1),
        FN_pct     = round(fn / n_total * 100, 1),
        stringsAsFactors = FALSE
      )
    })
  })

  do.call(rbind, unlist(result_list, recursive = FALSE))
}


# =============================================================================
# DISPERSION METRICS
# =============================================================================

#' Compute dispersion metrics by Comparison x Species
#'
#' Calculates MED, SD, CV, MAD, RCV, IQR, trimmed SD/CV for significant
#' proteins with correct direction.
#'
#' @param de_res Data frame with DE results (must include Species)
#' @param ev Data frame with expected values
#' @param alpha Significance threshold (default: 0.05)
#' @param lfc_thr Log fold-change threshold (default: 0)
#' @param p_col P-value column name (default: "adj.P.Val")
#' @param comparisons Comparisons to include (NULL = all)
#' @param assay Assay to filter (NULL = all)
#' @param trim Trimming proportion for trimmed SD/CV (default: 0.1)
#'
#' @return Data frame with dispersion metrics
#' @export
compute_dispersion_metrics <- function(de_res, ev,
                                       alpha = 0.05,
                                       lfc_thr = 0,
                                       p_col = "adj.P.Val",
                                       comparisons = NULL,
                                       assay = NULL,
                                       trim = 0.1) {
  # Prepare data
  prep <- .prepare_benchmark_data(de_res, ev, alpha, lfc_thr, p_col,
                                  comparisons, assay)
  de_res <- prep$de_res
  ev <- prep$expected_values
  p_col <- prep$p_col

  # Get expected species info
  comps <- unique(ev$Comparison)
  expected_lfc_map <- setNames(
    paste(ev$Species, ev$Comparison, sep = "||"),
    ev$expected_logFC
  )

  result_list <- lapply(comps, function(comp) {
    ev_comp <- ev[ev$Comparison == comp, , drop = FALSE]
    de_comp <- de_res[de_res$Comparison == comp, , drop = FALSE]
    all_species <- unique(de_comp$Species)

    lapply(all_species, function(sp) {
      de_sp <- de_comp[de_comp$Species == sp, , drop = FALSE]
      if (nrow(de_sp) == 0) return(NULL)

      # Get expected logFC for this species (NA if negative species)
      ev_row <- ev_comp[ev_comp$Species == sp, , drop = FALSE]
      exp_lfc <- if (nrow(ev_row) > 0) ev_row$expected_logFC[1] else NA_real_
      is_positive <- !is.na(exp_lfc)

      # Filter significant proteins
      pvals <- suppressWarnings(as.numeric(de_sp[[p_col]]))
      lfc <- as.numeric(de_sp$logFC)
      is_sig <- !is.na(pvals) & pvals <= alpha & abs(lfc) >= lfc_thr

      if (is_positive) {
        # Filter to correct direction
        correct_dir <- sign(lfc) == sign(exp_lfc)
        sig_correct <- is_sig & correct_dir
        lfc_values <- lfc[sig_correct]
      } else {
        # For negative species, all significant proteins
        lfc_values <- lfc[is_sig]
      }

      n_sig <- length(lfc_values)

      if (n_sig == 0) {
        return(data.frame(
          Comparison   = comp,
          Species      = sp,
          expected_logFC = ifelse(is_positive, exp_lfc, NA_real_),
          N_significant = 0L,
          MED = NA_real_, SD = NA_real_, CV = NA_real_,
          MAD = NA_real_, RCV = NA_real_, IQR = NA_real_,
          Trimmed_SD = NA_real_, Trimmed_CV = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      med_val <- median(lfc_values, na.rm = TRUE)
      sd_val  <- sd(lfc_values, na.rm = TRUE)
      mean_val <- mean(lfc_values, na.rm = TRUE)
      cv_val  <- if (abs(mean_val) > .Machine$double.eps) abs(sd_val / mean_val) * 100 else NA_real_
      mad_val <- mad(lfc_values, na.rm = TRUE)
      rcv_val <- if (abs(med_val) > .Machine$double.eps) abs(mad_val / med_val) * 100 else NA_real_
      iqr_val <- IQR(lfc_values, na.rm = TRUE)

      tsd_cv <- .trimmed_sd_cv(lfc_values, trim)

      data.frame(
        Comparison     = comp,
        Species        = sp,
        expected_logFC = ifelse(is_positive, exp_lfc, NA_real_),
        N_significant  = n_sig,
        MED = round(med_val, 4),
        SD  = round(sd_val, 4),
        CV  = round(cv_val, 2),
        MAD = round(mad_val, 4),
        RCV = round(rcv_val, 2),
        IQR = round(iqr_val, 4),
        Trimmed_SD = round(tsd_cv$trimmed_sd, 4),
        Trimmed_CV = round(tsd_cv$trimmed_cv, 2),
        stringsAsFactors = FALSE
      )
    })
  })

  do.call(rbind, unlist(result_list, recursive = FALSE))
}


# =============================================================================
# GGPLOT2 VISUALIZATIONS
# =============================================================================

#' Performance Heatmap (ggplot2)
#'
#' Heatmap of Comparisons x Metrics (AUC, Sensitivity, Specificity,
#' Precision, F1, Accuracy). Gradient red -> yellow -> green.
#'
#' @param metrics_table Data frame from compute_benchmark_metrics()
#' @param metrics_to_show Character vector of metrics to display
#' @param title Plot title
#' @param text_size Size of cell text labels
#' @param axis_text_size Size of axis text
#'
#' @return ggplot2 object
#' @export
benchmark_heatmap_gg <- function(
    metrics_table,
    metrics_to_show = c("AUC", "Sensitivity", "Specificity",
                        "Precision", "F1", "Accuracy"),
    title = "Benchmark Performance Heatmap",
    text_size = 4,
    axis_text_size = 11
) {
  # Filter to available metrics
  available <- intersect(metrics_to_show, names(metrics_table))
  if (length(available) == 0) {
    stop("Ninguna de las metricas solicitadas esta disponible en metrics_table")
  }

  # Reshape to long format
  plot_data <- metrics_table[, c("Comparison", available), drop = FALSE]

  plot_long <- tidyr::pivot_longer(
    plot_data,
    cols = -Comparison,
    names_to = "Metric",
    values_to = "Value"
  )

  # Order metrics as specified
  plot_long$Metric <- factor(plot_long$Metric, levels = rev(available))
  plot_long$Comparison <- factor(plot_long$Comparison,
                                 levels = unique(metrics_table$Comparison))

  # Build plot
  gg <- ggplot2::ggplot(plot_long,
                        ggplot2::aes(x = Comparison, y = Metric, fill = Value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(Value), "NA", sprintf("%.3f", Value))),
      size = text_size, color = "black", fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#E63946", mid = "#F4D35E", high = "#1a9850",
      midpoint = 0.5, limits = c(0, 1),
      na.value = "#CCCCCC",
      name = "Value"
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 15, color = "#1D3557"
      ),
      axis.text.x = ggplot2::element_text(
        size = axis_text_size, angle = 45, hjust = 1, color = "#495057"
      ),
      axis.text.y = ggplot2::element_text(
        size = axis_text_size, color = "#495057"
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  gg
}

#' Confusion Matrix Heatmap (ggplot2)
#'
#' Heatmap of Comparisons x (TP%, FP%, FN%, TN%) with white-to-blue gradient.
#'
#' @param confusion_df Data frame from .confusion_by_species()
#' @param title Plot title
#' @param text_size Size of cell text labels
#' @param axis_text_size Size of axis text
#'
#' @return ggplot2 object
#' @export
benchmark_confusion_gg <- function(
    confusion_df,
    title = "Confusion Matrix by Comparison x Species",
    text_size = 3.5,
    axis_text_size = 10
) {
  # Create label combining Comparison + Species
  confusion_df$Label <- paste0(confusion_df$Comparison, "\n", confusion_df$Species)

  # Reshape percentages to long format
  pct_cols <- c("TP_pct", "FP_pct", "FN_pct", "TN_pct")
  available_pct <- intersect(pct_cols, names(confusion_df))

  plot_long <- tidyr::pivot_longer(
    confusion_df[, c("Label", available_pct), drop = FALSE],
    cols = -Label,
    names_to = "Category",
    values_to = "Percentage"
  )

  # Clean category names
  plot_long$Category <- gsub("_pct$", "", plot_long$Category)
  plot_long$Category <- factor(plot_long$Category, levels = c("TP", "FP", "FN", "TN"))

  # Color mapping per category
  fill_colors <- c(
    "TP" = "#2A9D8F",
    "FP" = "#E63946",
    "FN" = "#F4A261",
    "TN" = "#457B9D"
  )

  gg <- ggplot2::ggplot(plot_long,
                        ggplot2::aes(x = Category, y = Label, fill = Category)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", Percentage)),
      size = text_size, color = "white", fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = fill_colors, name = "Classification") +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 14, color = "#1D3557"
      ),
      axis.text.x = ggplot2::element_text(
        size = axis_text_size, face = "bold", color = "#495057"
      ),
      axis.text.y = ggplot2::element_text(
        size = axis_text_size, color = "#495057"
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  gg
}

#' AUC Bar Chart (ggplot2)
#'
#' Horizontal bar chart of AUC per comparison, colored by value,
#' with a reference line at 0.5.
#'
#' @param metrics_table Data frame from compute_benchmark_metrics()
#' @param title Plot title
#' @param bar_width Bar width (default: 0.7)
#'
#' @return ggplot2 object
#' @export
benchmark_auc_bars_gg <- function(
    metrics_table,
    title = "AUC by Comparison",
    bar_width = 0.7
) {
  if (!"AUC" %in% names(metrics_table)) {
    stop("Columna 'AUC' no encontrada en metrics_table")
  }

  plot_data <- metrics_table[, c("Comparison", "AUC"), drop = FALSE]
  plot_data <- plot_data[!is.na(plot_data$AUC), , drop = FALSE]

  if (nrow(plot_data) == 0) {
    warning("No hay valores de AUC disponibles (requiere pROC)")
    return(NULL)
  }

  plot_data$Comparison <- factor(plot_data$Comparison,
                                 levels = rev(plot_data$Comparison))

  gg <- ggplot2::ggplot(plot_data,
                        ggplot2::aes(x = Comparison, y = AUC, fill = AUC)) +
    ggplot2::geom_col(width = bar_width) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", AUC)),
      hjust = -0.1, size = 3.5, fontface = "bold", color = "#1D3557"
    ) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                        color = "#E63946", linewidth = 0.6) +
    ggplot2::scale_fill_gradient2(
      low = "#E63946", mid = "#F4D35E", high = "#1a9850",
      midpoint = 0.75, limits = c(0, 1),
      name = "AUC"
    ) +
    ggplot2::coord_flip(ylim = c(0, max(plot_data$AUC, na.rm = TRUE) * 1.15)) +
    ggplot2::labs(
      title = title,
      subtitle = "Dashed line: Random classifier (AUC = 0.5)",
      x = NULL,
      y = "Area Under the Curve (AUC)"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 15, color = "#1D3557"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, size = 10, color = "#E63946", face = "italic"
      ),
      axis.text = ggplot2::element_text(size = 11, color = "#495057"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  gg
}

#' Grouped Metrics Bar Chart (ggplot2)
#'
#' Grouped bar chart: Comparisons x multiple metrics side by side.
#'
#' @param metrics_table Data frame from compute_benchmark_metrics()
#' @param metrics_to_show Character vector of metrics to display
#' @param title Plot title
#' @param bar_width Bar width (default: 0.7)
#'
#' @return ggplot2 object
#' @export
benchmark_metrics_bars_gg <- function(
    metrics_table,
    metrics_to_show = c("Sensitivity", "Specificity", "Precision", "F1", "Accuracy"),
    title = "Benchmark Metrics by Comparison",
    bar_width = 0.7
) {
  available <- intersect(metrics_to_show, names(metrics_table))
  if (length(available) == 0) {
    stop("Ninguna de las metricas solicitadas esta disponible")
  }

  plot_data <- metrics_table[, c("Comparison", available), drop = FALSE]

  plot_long <- tidyr::pivot_longer(
    plot_data,
    cols = -Comparison,
    names_to = "Metric",
    values_to = "Value"
  )

  plot_long$Metric <- factor(plot_long$Metric, levels = available)
  plot_long$Comparison <- factor(plot_long$Comparison,
                                 levels = unique(metrics_table$Comparison))

  # Color palette for metrics
  metric_colors <- c(
    "Sensitivity" = "#2A9D8F",
    "Specificity" = "#457B9D",
    "Precision"   = "#E9C46A",
    "F1"          = "#E63946",
    "Accuracy"    = "#264653",
    "NPV"         = "#F4A261",
    "MCC"         = "#9B5DE5",
    "AUC"         = "#00BBF9"
  )

  used_colors <- metric_colors[available]
  # Fill any missing colors
  if (any(is.na(used_colors))) {
    fallback <- c("#56B4E9", "#F0E442", "#CC79A7", "#882255")
    na_idx <- which(is.na(used_colors))
    used_colors[na_idx] <- fallback[seq_along(na_idx)]
  }

  gg <- ggplot2::ggplot(plot_long,
                        ggplot2::aes(x = Comparison, y = Value, fill = Metric)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = bar_width) +
    ggplot2::scale_fill_manual(values = used_colors) +
    ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Value",
      fill = "Metric"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 15, color = "#1D3557"
      ),
      axis.text.x = ggplot2::element_text(
        size = 11, angle = 45, hjust = 1, color = "#495057"
      ),
      axis.text.y = ggplot2::element_text(size = 11, color = "#495057"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )

  gg
}

#' Significant Proteins Stacked Bars (ggplot2)
#'
#' Stacked bar chart of significant proteins by species,
#' faceted by direction (UP/DOWN).
#'
#' @param classified_df Classified data frame
#' @param ev Expected values data frame
#' @param species_colors Named vector of colors per species (optional)
#' @param title Plot title
#'
#' @return ggplot2 object
#' @export
benchmark_signif_bars_gg <- function(
    classified_df,
    ev,
    species_colors = NULL,
    title = "Significant Proteins by Species and Direction"
) {
  # Filter significant proteins
  sig_df <- classified_df[classified_df$predicted == 1, , drop = FALSE]

  if (nrow(sig_df) == 0) {
    warning("No hay proteinas significativas para graficar")
    return(NULL)
  }

  # Assign direction based on logFC
  sig_df$Direction <- ifelse(sig_df$logFC > 0, "UP", "DOWN")

  # Count by Comparison x Species x Direction
  count_df <- as.data.frame(
    table(
      Comparison = sig_df$Comparison,
      Species    = sig_df$Species,
      Direction  = sig_df$Direction
    ),
    stringsAsFactors = FALSE
  )
  names(count_df)[4] <- "Count"
  count_df <- count_df[count_df$Count > 0, , drop = FALSE]

  if (nrow(count_df) == 0) {
    warning("No hay datos para graficar tras el conteo")
    return(NULL)
  }

  # Configure species colors
  all_species <- unique(count_df$Species)
  sp_colors <- .configure_species_colors(all_species, species_colors)

  gg <- ggplot2::ggplot(count_df,
                        ggplot2::aes(x = Comparison, y = Count, fill = Species)) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::facet_wrap(~ Direction) +
    ggplot2::scale_fill_manual(values = sp_colors) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Number of Proteins",
      fill = "Species"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 15, color = "#1D3557"
      ),
      axis.text.x = ggplot2::element_text(
        size = 11, angle = 45, hjust = 1, color = "#495057"
      ),
      axis.text.y = ggplot2::element_text(size = 11, color = "#495057"),
      strip.text = ggplot2::element_text(
        face = "bold", size = 12, color = "#1D3557"
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )

  gg
}


# =============================================================================
# ROC CURVES (ggplot2 via pROC)
# =============================================================================

#' ROC Curves by Comparison (ggplot2)
#'
#' Generates ROC curves for each comparison using pROC::ggroc().
#' Requires the pROC package.
#'
#' @param classified_df Classified data frame (output of .classify_all_comparisons)
#' @param p_col P-value column name used as predictor score
#' @param comparisons Comparisons to include (NULL = all)
#' @param title Plot title
#' @param zoom If TRUE, zoom into the low-FPR region (x = 0 to 0.1)
#' @param palette RColorBrewer palette name (default: "Set1")
#'
#' @return ggplot2 object or NULL if pROC is not available
#' @export
benchmark_roc_gg <- function(
    classified_df,
    p_col = "adj.P.Val",
    comparisons = NULL,
    title = "ROC Curves by Comparison",
    zoom = FALSE,
    palette = "Set1"
) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    warning("Paquete 'pROC' no instalado. No se pueden generar curvas ROC.\n",
            "Instalar con: install.packages('pROC')")
    return(NULL)
  }

  # Compute score: -log10(p-value) — higher = more significant
  p_col <- .detect_p_col(classified_df, p_col)
  classified_df$score <- -log10(pmax(as.numeric(classified_df[[p_col]]), 1e-300))

  # Filter comparisons
  if (!is.null(comparisons)) {
    classified_df <- classified_df[classified_df$Comparison %in% comparisons, , drop = FALSE]
  }

  # Filter valid rows
  classified_df <- classified_df[!is.na(classified_df$truth) & !is.na(classified_df$score), ,
                                 drop = FALSE]

  if (nrow(classified_df) == 0) {
    warning("No hay datos validos para generar curvas ROC")
    return(NULL)
  }

  # Build ROC objects per comparison
  comp_list <- split(classified_df, classified_df$Comparison)

  roc_list <- lapply(comp_list, function(d) {
    if (length(unique(d$truth)) < 2 || nrow(d) < 10) return(NULL)
    tryCatch(
      pROC::roc(response = d$truth, predictor = d$score, quiet = TRUE),
      error = function(e) NULL
    )
  })

  roc_list <- Filter(Negate(is.null), roc_list)

  if (length(roc_list) == 0) {
    warning("No se pudieron generar curvas ROC para ninguna comparacion")
    return(NULL)
  }

  # AUC labels for legend
  aucs <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))
  labels <- paste0(names(roc_list), " (AUC=", sprintf("%.3f", aucs), ")")

  # Build plot
  gg <- pROC::ggroc(roc_list, legacy.axes = TRUE, linewidth = 1) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "gray50", linewidth = 0.5
    ) +
    ggplot2::labs(
      title = title,
      subtitle = if (zoom) {
        "Zoom: low False Positive Rate region (0-10%)"
      } else {
        "Diagonal = random classifier"
      },
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Comparison"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 15, color = "#1D3557"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, size = 11, color = "#6C757D"
      ),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Apply palette
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    n_colors <- min(length(roc_list), RColorBrewer::brewer.pal.info[palette, "maxcolors"])
    pal_colors <- RColorBrewer::brewer.pal(max(n_colors, 3), palette)
    gg <- gg + ggplot2::scale_color_manual(values = pal_colors, labels = labels)
  } else {
    gg <- gg + ggplot2::scale_color_discrete(labels = labels)
  }

  # Zoom or full view
  if (zoom) {
    gg <- gg + ggplot2::coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 1))
  } else {
    gg <- gg + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  }

  gg
}


# =============================================================================
# BENCHMARK VOLCANO PLOT — HIGHCHARTER
# =============================================================================

#' Benchmark Volcano Plot (Highcharter)
#'
#' Volcano plot colored by Species with plotLines at expected logFC values,
#' dispersion stats in subtitle.
#'
#' @param de_res_comp DE results for one comparison
#' @param ev_comp Expected values for one comparison
#' @param disp_comp Dispersion metrics for one comparison (optional)
#' @param alpha Significance threshold
#' @param lfc_thr Log fold-change threshold
#' @param p_col P-value column
#' @param species_colors Named vector of colors per species
#' @param point_size Marker radius (default: 4)
#' @param title Custom title (optional)
#' @param height Chart height in pixels (optional)
#'
#' @return Highchart object
#' @export
benchmark_volcano_hc <- function(
    de_res_comp,
    ev_comp,
    disp_comp = NULL,
    alpha = 0.05,
    lfc_thr = 0,
    p_col = "adj.P.Val",
    species_colors = NULL,
    point_size = 4,
    title = NULL,
    height = NULL
) {
  if (!requireNamespace("highcharter", quietly = TRUE)) {
    stop("Se requiere el paquete 'highcharter'")
  }

  p_col <- .detect_p_col(de_res_comp, p_col)

  # Compute -log10(p)
  de_res_comp$minusLog10P <- -log10(as.numeric(de_res_comp[[p_col]]))
  max_finite <- max(de_res_comp$minusLog10P[is.finite(de_res_comp$minusLog10P)],
                    na.rm = TRUE)
  de_res_comp$minusLog10P[!is.finite(de_res_comp$minusLog10P)] <- max_finite * 1.1

  de_res_comp$pval_fmt <- sprintf("%.3g", as.numeric(de_res_comp[[p_col]]))

  # Significance
  pvals <- suppressWarnings(as.numeric(de_res_comp[[p_col]]))
  lfc <- as.numeric(de_res_comp$logFC)
  is_sig <- !is.na(pvals) & pvals <= alpha & abs(lfc) >= lfc_thr

  # All species
  all_species <- unique(de_res_comp$Species)
  sp_colors <- .configure_species_colors(all_species, species_colors)

  # Comparison label
  comp <- unique(de_res_comp$Comparison)[1]

  # Title
  if (is.null(title)) {
    assay_label <- ""
    if ("Assay" %in% names(de_res_comp)) {
      assays <- unique(de_res_comp$Assay)
      if (length(assays) == 1) assay_label <- paste0(" (", assays, ")")
    }
    title <- paste0("Benchmark Volcano: ", comp, assay_label)
  }

  # --- Build subtitle with dispersion stats ---
  subtitle_html <- ""
  if (!is.null(disp_comp) && nrow(disp_comp) > 0) {
    subtitle_parts <- vapply(seq_len(nrow(disp_comp)), function(i) {
      sp <- disp_comp$Species[i]
      sp_col <- sp_colors[sp] %||% "#666666"
      med <- disp_comp$MED[i]
      mad <- disp_comp$MAD[i]
      rcv <- disp_comp$RCV[i]
      n <- disp_comp$N_significant[i]

      if (is.na(med)) return("")

      sprintf(
        "<span style='color:%s;font-weight:bold;'>%s</span>: MED=%.2f MAD=%.3f RCV=%.1f%% (n=%d)",
        sp_col, sp, med, mad, rcv, n
      )
    }, character(1))

    subtitle_parts <- subtitle_parts[subtitle_parts != ""]
    subtitle_html <- paste(subtitle_parts, collapse = " | ")
  }

  # --- Build series per species group ---
  # Groups: {Species}_Significant (color) + Non_Significant (grey)
  positive_species <- ev_comp$Species

  # Create series list
  series_list <- list()

  for (sp in all_species) {
    sp_mask <- de_res_comp$Species == sp
    sp_data <- de_res_comp[sp_mask, , drop = FALSE]
    sp_sig <- sp_data[is_sig[sp_mask], , drop = FALSE]
    sp_ns <- sp_data[!is_sig[sp_mask], , drop = FALSE]

    sp_color <- unname(sp_colors[sp])

    # Gene names column
    gene_col <- if ("Gene.Names" %in% names(sp_data)) "Gene.Names" else "Protein.IDs"

    # Significant series
    if (nrow(sp_sig) > 0) {
      sig_points <- lapply(seq_len(nrow(sp_sig)), function(i) {
        list(
          x = sp_sig$logFC[i],
          y = sp_sig$minusLog10P[i],
          gene = as.character(sp_sig[[gene_col]][i]),
          protein = as.character(sp_sig$Protein.IDs[i]),
          pval = sp_sig$pval_fmt[i],
          species = sp
        )
      })

      series_list <- c(series_list, list(list(
        name = paste0(sp, " (Significant)"),
        type = "scatter",
        data = sig_points,
        color = sp_color,
        marker = list(
          radius = point_size,
          symbol = "circle",
          lineWidth = 0
        ),
        showInLegend = TRUE,
        enableMouseTracking = TRUE
      )))
    }

    # Non-significant series
    if (nrow(sp_ns) > 0) {
      ns_points <- lapply(seq_len(nrow(sp_ns)), function(i) {
        list(
          x = sp_ns$logFC[i],
          y = sp_ns$minusLog10P[i],
          gene = as.character(sp_ns[[gene_col]][i]),
          protein = as.character(sp_ns$Protein.IDs[i]),
          pval = sp_ns$pval_fmt[i],
          species = sp
        )
      })

      ns_color <- .hex_to_rgba(sp_color, 0.25)

      series_list <- c(series_list, list(list(
        name = paste0(sp, " (NS)"),
        type = "scatter",
        data = ns_points,
        color = ns_color,
        marker = list(
          radius = point_size - 1,
          symbol = "circle",
          lineWidth = 0
        ),
        showInLegend = TRUE,
        enableMouseTracking = TRUE
      )))
    }
  }

  # --- xAxis plotLines: expected logFC ---
  x_plotlines <- list(
    list(value = 0, color = "#1D3557", width = 1, zIndex = 3)
  )

  if (lfc_thr > 0) {
    x_plotlines <- c(x_plotlines, list(
      list(value = -lfc_thr, color = "#6C757D", width = 1,
           dashStyle = "Dash", zIndex = 2),
      list(value = lfc_thr, color = "#6C757D", width = 1,
           dashStyle = "Dash", zIndex = 2)
    ))
  }

  # Expected logFC lines per species
  for (i in seq_len(nrow(ev_comp))) {
    sp <- ev_comp$Species[i]
    exp_lfc <- ev_comp$expected_logFC[i]
    sp_col <- unname(sp_colors[sp] %||% "#666666")

    x_plotlines <- c(x_plotlines, list(list(
      value = exp_lfc,
      color = sp_col,
      width = 2,
      dashStyle = "LongDash",
      zIndex = 4,
      label = list(
        text = sprintf("%s (%.2f)", sp, exp_lfc),
        style = list(
          color = sp_col,
          fontSize = "10px",
          fontWeight = "bold"
        ),
        align = "left",
        x = 5,
        y = 12
      )
    )))
  }

  # --- Build highchart ---
  hc <- highcharter::highchart() |>
    highcharter::hc_chart(
      type = "scatter",
      zoomType = "xy",
      backgroundColor = "#FFFFFF",
      style = list(fontFamily = "Inter, -apple-system, sans-serif"),
      height = height
    ) |>
    highcharter::hc_title(
      text = title,
      style = list(fontSize = "16px", fontWeight = "600", color = "#1D3557")
    ) |>
    highcharter::hc_xAxis(
      title = list(
        text = "log<sub>2</sub> Fold Change",
        useHTML = TRUE,
        style = list(fontSize = "13px", color = "#495057")
      ),
      gridLineWidth = 0,
      lineColor = "#DEE2E6",
      tickColor = "#DEE2E6",
      plotLines = x_plotlines
    ) |>
    highcharter::hc_yAxis(
      title = list(
        text = paste0("-log<sub>10</sub>(", p_col, ")"),
        useHTML = TRUE,
        style = list(fontSize = "13px", color = "#495057")
      ),
      lineColor = "#DEE2E6",
      lineWidth = 1,
      tickColor = "#DEE2E6",
      gridLineColor = "#F1F3F4",
      gridLineDashStyle = "Dot",
      plotLines = list(
        list(
          value = -log10(alpha),
          color = "#6C757D",
          width = 1,
          dashStyle = "Dash",
          zIndex = 2,
          label = list(
            text = paste0("\u03b1 = ", alpha),
            style = list(color = "#6C757D", fontSize = "10px"),
            align = "right",
            x = -10,
            y = 12
          )
        )
      )
    ) |>
    highcharter::hc_tooltip(
      useHTML = TRUE,
      backgroundColor = "rgba(255, 255, 255, 0.95)",
      borderColor = "#DEE2E6",
      borderRadius = 8,
      shadow = TRUE,
      style = list(fontSize = "12px"),
      headerFormat = "",
      pointFormat = paste0(
        "<div style='padding: 4px;'>",
        "<b style='font-size: 13px; color: #1D3557;'>{point.gene}</b><br/>",
        "<span style='color: #6C757D;'>Protein:</span> {point.protein}<br/>",
        "<span style='color: #6C757D;'>Species:</span> {point.species}<br/>",
        "<span style='color: #6C757D;'>log\u2082FC:</span> <b>{point.x:.3f}</b><br/>",
        "<span style='color: #6C757D;'>", p_col, ":</span> <b>{point.pval}</b>",
        "</div>"
      )
    ) |>
    highcharter::hc_legend(
      title = list(text = "", style = list(fontStyle = "normal")),
      layout = "horizontal",
      align = "center",
      verticalAlign = "bottom",
      itemStyle = list(fontSize = "11px", fontWeight = "normal")
    ) |>
    highcharter::hc_plotOptions(
      scatter = list(
        marker = list(
          states = list(
            hover = list(
              radiusPlus = 2,
              lineWidthPlus = 1,
              lineColor = "#1D3557"
            )
          )
        ),
        turboThreshold = nrow(de_res_comp) + 100
      )
    ) |>
    highcharter::hc_exporting(
      enabled = TRUE,
      buttons = list(
        contextButton = list(
          menuItems = c("downloadPNG", "downloadSVG", "downloadPDF")
        )
      )
    )

  # Add subtitle
  if (nchar(subtitle_html) > 0) {
    hc <- hc |>
      highcharter::hc_subtitle(
        text = subtitle_html,
        useHTML = TRUE,
        style = list(fontSize = "11px", color = "#6C757D")
      )
  }

  # Add data series
  for (series in series_list) {
    hc <- hc |>
      highcharter::hc_add_series(
        name = series$name,
        type = series$type,
        data = series$data,
        color = series$color,
        marker = series$marker,
        showInLegend = series$showInLegend,
        enableMouseTracking = series$enableMouseTracking
      )
  }

  hc
}

#' List of Benchmark Volcano Plots (Highcharter)
#'
#' Generates one volcano plot per comparison.
#'
#' @param de_res Data frame with DE results
#' @param ev Data frame with expected values
#' @param disp Dispersion metrics data frame (optional)
#' @param alpha Significance threshold (default: 0.05)
#' @param lfc_thr Log fold-change threshold (default: 0)
#' @param p_col P-value column (default: "adj.P.Val")
#' @param comparisons Comparisons to include (NULL = all)
#' @param assay Assay to filter (NULL = all)
#' @param species_colors Named vector of colors per species (optional)
#' @param point_size Marker radius (default: 4)
#' @param height Chart height in pixels (optional)
#'
#' @return Named list of highchart objects
#' @export
benchmark_volcano_hc_list <- function(
    de_res,
    ev,
    disp = NULL,
    alpha = 0.05,
    lfc_thr = 0,
    p_col = "adj.P.Val",
    comparisons = NULL,
    assay = NULL,
    species_colors = NULL,
    point_size = 4,
    height = NULL
) {
  # Prepare data
  prep <- .prepare_benchmark_data(de_res, ev, alpha, lfc_thr, p_col,
                                  comparisons, assay)
  de_res <- prep$de_res
  ev <- prep$expected_values
  p_col <- prep$p_col
  comps <- prep$comparisons

  # Configure species colors once
  all_species <- unique(de_res$Species)
  sp_colors <- .configure_species_colors(all_species, species_colors)

  hc_list <- lapply(comps, function(comp) {
    de_comp <- de_res[de_res$Comparison == comp, , drop = FALSE]
    ev_comp <- ev[ev$Comparison == comp, , drop = FALSE]

    disp_comp <- NULL
    if (!is.null(disp)) {
      disp_comp <- disp[disp$Comparison == comp, , drop = FALSE]
    }

    if (nrow(de_comp) == 0 || nrow(ev_comp) == 0) return(NULL)

    tryCatch({
      benchmark_volcano_hc(
        de_res_comp = de_comp,
        ev_comp = ev_comp,
        disp_comp = disp_comp,
        alpha = alpha,
        lfc_thr = lfc_thr,
        p_col = p_col,
        species_colors = sp_colors,
        point_size = point_size,
        height = height
      )
    }, error = function(e) {
      warning(sprintf("Error generando volcano para %s: %s", comp, e$message))
      NULL
    })
  })

  names(hc_list) <- comps
  Filter(Negate(is.null), hc_list)
}


# =============================================================================
# EXPORT
# =============================================================================

#' Export benchmark data to TSV or Parquet
#'
#' @param data Data frame to export
#' @param filepath Full file path (with extension)
#' @param format "tsv" or "parquet"
#' @keywords internal
.export_benchmark_data <- function(data, filepath, format = "tsv") {
  if (format == "tsv") {
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_tsv(data, filepath)
    } else {
      write.table(data, filepath, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  } else if (format == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      warning("Paquete 'arrow' no instalado. Exportando como TSV.")
      filepath <- sub("\\.parquet$", ".tsv", filepath)
      write.table(data, filepath, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      arrow::write_parquet(data, filepath)
    }
  }

  invisible(filepath)
}

#' Export ggplot2 chart to PNG
#'
#' @param gg ggplot2 object
#' @param filepath Full file path
#' @param width Width in inches (default: 10)
#' @param height Height in inches (default: 7)
#' @param dpi Resolution (default: 300)
#' @keywords internal
.export_gg_plot <- function(gg, filepath, width = 10, height = 7, dpi = 300) {
  if (is.null(gg)) return(invisible(NULL))

  tryCatch({
    ggplot2::ggsave(filepath, plot = gg, width = width, height = height,
                    dpi = dpi, bg = "white")
  }, error = function(e) {
    warning("Error exportando grafico a ", filepath, ": ", e$message)
  })

  invisible(filepath)
}


# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Benchmark Proteomics Differential Expression Results
#'
#' Evaluates the reliability of DE results by comparing them against
#' known ground truth (multi-species spike-in experiments). Computes
#' classification metrics, dispersion statistics, and generates
#' visualizations.
#'
#' @param de_res Data frame with DE results. Required columns:
#'   Protein.IDs, logFC, P.Value, adj.P.Val, Change, Comparison, Species.
#'   Optional: Gene.Names, Assay.
#' @param expected_values Data frame with expected values. Required columns:
#'   Comparison, Species, expected_logFC.
#' @param alpha Significance threshold (default: 0.05)
#' @param lfc_thr Log fold-change threshold (default: 0)
#' @param p_col P-value column to use (default: "adj.P.Val")
#' @param comparisons Comparisons to include (NULL = all)
#' @param assay Assay to filter (NULL = all)
#' @param output_dir Directory for exported files (NULL = no export)
#' @param verbose Print progress messages (default: TRUE)
#' @param species_colors Named vector of colors per species (optional)
#'
#' @return List with:
#'   \itemize{
#'     \item metrics_table: Classification metrics per comparison
#'     \item confusion_by_species: Confusion matrix per Comparison x Species
#'     \item dispersion_metrics: Dispersion stats per Comparison x Species
#'     \item classified_df: Full classified data frame
#'     \item gg_heatmap: ggplot2 performance heatmap
#'     \item gg_confusion: ggplot2 confusion matrix heatmap
#'     \item gg_auc_bars: ggplot2 AUC bar chart
#'     \item gg_metrics_bars: ggplot2 grouped metrics bar chart
#'     \item gg_signif_bars: ggplot2 significant proteins stacked bars
#'     \item hc_volcano_list: Named list of Highcharter volcano plots
#'     \item parameters: List of parameters used
#'   }
#'
#' @examples
#' \dontrun{
#' expected <- data.frame(
#'   Comparison   = c("B/A", "B/A", "C/A", "C/A"),
#'   Species      = c("ECOLI", "YEAST", "ECOLI", "YEAST"),
#'   expected_logFC = c(1.0, -0.58, 1.58, -1.60)
#' )
#'
#' result <- benchmarking_proteomics(
#'   de_res = DEPs_results_with_species,
#'   expected_values = expected,
#'   alpha = 0.05,
#'   output_dir = "data/benchmark"
#' )
#'
#' result$gg_heatmap
#' result$hc_volcano_list[["B/A"]]
#' }
#'
#' @export
benchmarking_proteomics <- function(
    de_res,
    expected_values,
    alpha = 0.05,
    lfc_thr = 0,
    p_col = "adj.P.Val",
    comparisons = NULL,
    assay = NULL,
    output_dir = NULL,
    verbose = TRUE,
    species_colors = NULL
) {
  # === STEP 1: Prepare data ===
  if (verbose) cat("\n=== BENCHMARKING PROTEOMICS ===\n")

  prep <- .prepare_benchmark_data(de_res, expected_values, alpha, lfc_thr,
                                  p_col, comparisons, assay)
  de_res <- prep$de_res
  ev <- prep$expected_values
  p_col <- prep$p_col
  comps <- prep$comparisons

  if (verbose) {
    cat("- Comparaciones:", paste(comps, collapse = ", "), "\n")
    cat("- Especies en datos:", paste(unique(de_res$Species), collapse = ", "), "\n")
    cat("- Alpha:", alpha, "| LFC threshold:", lfc_thr, "\n")
    cat("- P-value column:", p_col, "\n")
  }

  # Configure species colors
  all_species <- unique(de_res$Species)
  sp_colors <- .configure_species_colors(all_species, species_colors)

  # === STEP 2: Classification ===
  if (verbose) cat("\n--- Clasificacion ground truth ---\n")
  classified_df <- .classify_all_comparisons(de_res, ev, alpha, lfc_thr, p_col)

  if (verbose) {
    for (comp in comps) {
      df_comp <- classified_df[classified_df$Comparison == comp, , drop = FALSE]
      tp <- sum(df_comp$classification == "TP")
      fp <- sum(df_comp$classification == "FP")
      tn <- sum(df_comp$classification == "TN")
      fn <- sum(df_comp$classification == "FN")
      cat(sprintf("  %s: TP=%d FP=%d TN=%d FN=%d\n", comp, tp, fp, tn, fn))
    }
  }

  # === STEP 3: Metrics ===
  if (verbose) cat("\n--- Metricas de clasificacion ---\n")
  metrics_table <- compute_benchmark_metrics(de_res, ev, alpha, lfc_thr, p_col,
                                             comparisons = comps, assay = assay)

  if (!requireNamespace("pROC", quietly = TRUE) && verbose) {
    cat("  [NOTA] Paquete 'pROC' no instalado. AUC = NA.\n")
    cat("  Instalar con: install.packages('pROC')\n")
  }

  if (verbose) {
    for (i in seq_len(nrow(metrics_table))) {
      cat(sprintf("  %s: Sens=%.3f Spec=%.3f Prec=%.3f F1=%.3f AUC=%s\n",
                  metrics_table$Comparison[i],
                  metrics_table$Sensitivity[i],
                  metrics_table$Specificity[i],
                  metrics_table$Precision[i],
                  metrics_table$F1[i],
                  ifelse(is.na(metrics_table$AUC[i]), "NA",
                         sprintf("%.3f", metrics_table$AUC[i]))))
    }
  }

  # === STEP 4: Confusion by species ===
  confusion_by_species_df <- .confusion_by_species(classified_df)

  # === STEP 5: Dispersion ===
  if (verbose) cat("\n--- Metricas de dispersion ---\n")
  dispersion_df <- compute_dispersion_metrics(de_res, ev, alpha, lfc_thr, p_col,
                                              comparisons = comps, assay = assay)

  if (verbose) {
    pos_disp <- dispersion_df[!is.na(dispersion_df$expected_logFC), , drop = FALSE]
    for (i in seq_len(nrow(pos_disp))) {
      cat(sprintf("  %s | %s: MED=%.3f MAD=%.4f RCV=%.1f%% (n=%d)\n",
                  pos_disp$Comparison[i],
                  pos_disp$Species[i],
                  pos_disp$MED[i],
                  pos_disp$MAD[i],
                  pos_disp$RCV[i],
                  pos_disp$N_significant[i]))
    }
  }

  # === STEP 6: Visualizations ===
  if (verbose) cat("\n--- Generando visualizaciones ---\n")

  # ggplot2 charts
  gg_heatmap <- tryCatch({
    if (verbose) cat("  - Heatmap de metricas (ggplot2)\n")
    benchmark_heatmap_gg(metrics_table)
  }, error = function(e) {
    warning("Error generando heatmap: ", e$message)
    NULL
  })

  gg_confusion <- tryCatch({
    if (verbose) cat("  - Heatmap de confusion (ggplot2)\n")
    benchmark_confusion_gg(confusion_by_species_df)
  }, error = function(e) {
    warning("Error generando confusion heatmap: ", e$message)
    NULL
  })

  gg_auc_bars <- tryCatch({
    if (verbose) cat("  - Barras AUC (ggplot2)\n")
    benchmark_auc_bars_gg(metrics_table)
  }, error = function(e) {
    warning("Error generando AUC bars: ", e$message)
    NULL
  })

  gg_metrics_bars <- tryCatch({
    if (verbose) cat("  - Barras de metricas agrupadas (ggplot2)\n")
    benchmark_metrics_bars_gg(metrics_table)
  }, error = function(e) {
    warning("Error generando metrics bars: ", e$message)
    NULL
  })

  gg_signif_bars <- tryCatch({
    if (verbose) cat("  - Barras de significativas por especie (ggplot2)\n")
    benchmark_signif_bars_gg(classified_df, ev, species_colors = sp_colors)
  }, error = function(e) {
    warning("Error generando signif bars: ", e$message)
    NULL
  })

  # Highcharter volcano
  hc_volcano_list <- tryCatch({
    if (verbose) cat("  - Volcano plots benchmark (Highcharter)\n")
    benchmark_volcano_hc_list(
      de_res, ev, disp = dispersion_df,
      alpha = alpha, lfc_thr = lfc_thr, p_col = p_col,
      comparisons = comps, assay = assay,
      species_colors = sp_colors, point_size = 4
    )
  }, error = function(e) {
    warning("Error generando volcano plots: ", e$message)
    list()
  })

  if (verbose) cat("  - Total volcanos generados:", length(hc_volcano_list), "\n")

  # ROC curves (require pROC)
  gg_roc <- tryCatch({
    if (verbose) cat("  - Curvas ROC (ggplot2 + pROC)\n")
    benchmark_roc_gg(classified_df, p_col = p_col, comparisons = comps)
  }, error = function(e) {
    warning("Error generando curvas ROC: ", e$message)
    NULL
  })

  gg_roc_zoom <- tryCatch({
    if (verbose) cat("  - Curvas ROC zoom (ggplot2 + pROC)\n")
    benchmark_roc_gg(classified_df, p_col = p_col, comparisons = comps,
                     title = "ROC Curves by Comparison (Zoom)",
                     zoom = TRUE)
  }, error = function(e) {
    warning("Error generando curvas ROC zoom: ", e$message)
    NULL
  })

  # === STEP 7: Performance summary ===
  metrics_table$Performance <- ifelse(
    is.na(metrics_table$F1), "NA",
    ifelse(metrics_table$F1 >= 0.9, "Excellent",
    ifelse(metrics_table$F1 >= 0.8, "Very Good",
    ifelse(metrics_table$F1 >= 0.7, "Good",
    ifelse(metrics_table$F1 >= 0.5, "Acceptable", "Poor"))))
  )

  # === STEP 8: Export ===
  if (!is.null(output_dir)) {
    if (verbose) cat("\n--- Exportando resultados ---\n")

    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      if (verbose) cat("  - Directorio creado:", output_dir, "\n")
    }

    # TSV exports
    .export_benchmark_data(
      metrics_table,
      file.path(output_dir, "benchmark_metrics.tsv"), "tsv"
    )
    if (verbose) cat("  - benchmark_metrics.tsv\n")

    .export_benchmark_data(
      confusion_by_species_df,
      file.path(output_dir, "benchmark_confusion_by_species.tsv"), "tsv"
    )
    if (verbose) cat("  - benchmark_confusion_by_species.tsv\n")

    .export_benchmark_data(
      dispersion_df,
      file.path(output_dir, "benchmark_dispersion.tsv"), "tsv"
    )
    if (verbose) cat("  - benchmark_dispersion.tsv\n")

    # Summary with Performance column
    summary_df <- metrics_table[, c("Comparison", "Sensitivity", "Specificity",
                                    "Precision", "F1", "AUC", "Accuracy",
                                    "MCC", "Performance"), drop = FALSE]
    .export_benchmark_data(
      summary_df,
      file.path(output_dir, "benchmark_summary.tsv"), "tsv"
    )
    if (verbose) cat("  - benchmark_summary.tsv\n")

    # PNG exports
    .export_gg_plot(gg_heatmap,
                    file.path(output_dir, "benchmark_heatmap.png"),
                    width = 10, height = 6)
    if (verbose) cat("  - benchmark_heatmap.png\n")

    .export_gg_plot(gg_confusion,
                    file.path(output_dir, "benchmark_confusion.png"),
                    width = 10, height = 7)
    if (verbose) cat("  - benchmark_confusion.png\n")

    .export_gg_plot(gg_auc_bars,
                    file.path(output_dir, "benchmark_auc_bars.png"),
                    width = 8, height = 5)
    if (verbose) cat("  - benchmark_auc_bars.png\n")

    .export_gg_plot(gg_metrics_bars,
                    file.path(output_dir, "benchmark_metrics_bars.png"),
                    width = 11, height = 6)
    if (verbose) cat("  - benchmark_metrics_bars.png\n")

    .export_gg_plot(gg_signif_bars,
                    file.path(output_dir, "benchmark_signif_bars.png"),
                    width = 10, height = 6)
    if (verbose) cat("  - benchmark_signif_bars.png\n")

    .export_gg_plot(gg_roc,
                    file.path(output_dir, "benchmark_roc.png"),
                    width = 10, height = 7)
    if (verbose) cat("  - benchmark_roc.png\n")

    .export_gg_plot(gg_roc_zoom,
                    file.path(output_dir, "benchmark_roc_zoom.png"),
                    width = 10, height = 7)
    if (verbose) cat("  - benchmark_roc_zoom.png\n")
  }

  if (verbose) cat("\n=== BENCHMARKING COMPLETADO ===\n\n")

  # === RETURN ===
  list(
    metrics_table        = metrics_table,
    confusion_by_species = confusion_by_species_df,
    dispersion_metrics   = dispersion_df,
    classified_df        = classified_df,
    gg_heatmap           = gg_heatmap,
    gg_confusion         = gg_confusion,
    gg_auc_bars          = gg_auc_bars,
    gg_metrics_bars      = gg_metrics_bars,
    gg_signif_bars       = gg_signif_bars,
    gg_roc               = gg_roc,
    gg_roc_zoom          = gg_roc_zoom,
    hc_volcano_list      = hc_volcano_list,
    parameters           = list(
      alpha          = alpha,
      lfc_thr        = lfc_thr,
      p_col          = p_col,
      comparisons    = comps,
      assay          = assay,
      species_colors = sp_colors
    )
  )
}


# =============================================================================
# EXAMPLES
# =============================================================================

# source("R/Benchmarking.R")
#
# # --- Define expected values (spike-in design) ---
# expected <- data.frame(
#   Comparison   = c("B/A", "B/A", "C/A", "C/A"),
#   Species      = c("ECOLI", "YEAST", "ECOLI", "YEAST"),
#   expected_logFC = c(1.0, -0.58, 1.58, -1.60)
# )
#
# # --- Run full benchmark ---
# result <- benchmarking_proteomics(
#   de_res          = DEPs_results_with_species,
#   expected_values = expected,
#   alpha           = 0.05,
#   lfc_thr         = 0,
#   output_dir      = "data/benchmark",
#   verbose         = TRUE
# )
#
# # --- Access results ---
# result$metrics_table
# result$dispersion_metrics
# result$confusion_by_species
#
# # --- View ggplot2 visualizations ---
# result$gg_heatmap
# result$gg_confusion
# result$gg_auc_bars
# result$gg_metrics_bars
# result$gg_signif_bars
# result$gg_roc
# result$gg_roc_zoom
#
# # --- View Highcharter volcanos ---
# result$hc_volcano_list[["B/A"]]
# result$hc_volcano_list[["C/A"]]
#
# # --- Use individual functions ---
# metrics <- compute_benchmark_metrics(
#   de_res = DEPs_results_with_species,
#   ev     = expected,
#   alpha  = 0.05
# )
#
# dispersion <- compute_dispersion_metrics(
#   de_res = DEPs_results_with_species,
#   ev     = expected,
#   alpha  = 0.05
# )
#
# # --- Custom heatmap ---
# benchmark_heatmap_gg(
#   metrics,
#   metrics_to_show = c("AUC", "F1", "MCC"),
#   title = "Custom Benchmark Heatmap"
# )
#
# # --- Custom volcano for a single comparison ---
# de_BA <- DEPs_results_with_species[
#   DEPs_results_with_species$Comparison == "B/A", ]
# ev_BA <- expected[expected$Comparison == "B/A", ]
# disp_BA <- dispersion[dispersion$Comparison == "B/A", ]
#
# benchmark_volcano_hc(
#   de_res_comp = de_BA,
#   ev_comp     = ev_BA,
#   disp_comp   = disp_BA,
#   alpha       = 0.05,
#   species_colors = c(ECOLI = "#D55E00", YEAST = "#0072B2", HUMAN = "#999999")
# )
