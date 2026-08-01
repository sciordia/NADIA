# =============================================================================
# tidyHeatmap heatmaps for proteomics data
# =============================================================================



# -----------------------------------------------------------------------------
# Wrapper and print method for heatmaps carrying a custom title
# -----------------------------------------------------------------------------

#' Wrap a Heatmap Together With Its Title
#'
#' @param hm tidyHeatmap object (InputHeatmap)
#' @param title Heatmap title
#' @param title_size Font size of the title
#' @param title_face Font face of the title
#'
#' @return proteomics_heatmap object (S3 list)
wrap_heatmap_with_title <- function(hm, title, title_size = 14, title_face = "bold") {

  structure(
    list(
      heatmap = hm,
      title = title,
      title_size = title_size,
      title_face = title_face
    ),
    class = "proteomics_heatmap"
  )
}

#' Draw a proteomics heatmap
#'
#' Renders the object built by [proteomics_heatmap()] on the current graphics
#' device, adding the title stored in the object. Printing is what actually draws
#' the heatmap: building it does not.
#'
#' @param x A `proteomics_heatmap` object, as returned by [proteomics_heatmap()].
#' @param ... Further arguments passed to `ComplexHeatmap::draw()`.
#' @return The drawn `HeatmapList`, invisibly. Called for its side effect.
#'
#' @examples
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("tidyHeatmap", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'   hm <- proteomics_heatmap(res$PCA_Input, mode = "any")
#'
#'   # Drawing needs a graphics device; send it to a temporary file
#'   f <- tempfile(fileext = ".png")
#'   grDevices::png(f)
#'   print(hm)
#'   grDevices::dev.off()
#'   file.remove(f)
#' }
#'
#' @export
print.proteomics_heatmap <- function(x, ...) {
  hm <- x$heatmap

  # Convert the InputHeatmap into a ComplexHeatmap using tidyHeatmap's method
  if (inherits(hm, "InputHeatmap")) {
    # Either use as.list to pull out the components and rebuild the object,
    # or simply convert it with tidyHeatmap's own internal method
    ht <- tryCatch({
      # Try the slot directly (older versions)
      methods::slot(hm, "ht")
    }, error = function(e) {
      tryCatch({
        # Try input_heatmap (newer versions)
        methods::slot(hm, "input_heatmap")
      }, error = function(e2) {
        # Fall back on tidyHeatmap's conversion to ComplexHeatmap
        tidyHeatmap::as_ComplexHeatmap(hm)
      })
    })
  } else {
    ht <- hm
  }

  # Draw with the title
  ComplexHeatmap::draw(
    ht,
    column_title = x$title,
    column_title_gp = grid::gpar(fontsize = x$title_size, fontface = x$title_face),
    ...
  )

  invisible(x)
}



# -----------------------------------------------------------------------------
# Helper that builds the colour palette for the heatmap
# -----------------------------------------------------------------------------

#' Build the Colour Palette for the Heatmap Values
#'
#' @param palette Palette specification:
#'   - NULL: use the default palette (RdBu)
#'   - Vector of colours: use those colours directly
#'   - String "paletteer::" (e.g. "viridis::viridis"): use paletteer
#'   - String "brewer:" (e.g. "brewer:RdBu"): use RColorBrewer
#' @param n Number of colours to generate
#' @param reverse Reverse the palette (default: FALSE)
#'
#' @return Vector of colours, or a colorRamp2 function
get_heatmap_palette <- function(palette = NULL,
                                n = 11,
                                reverse = FALSE) {

  # Default palette (diverging blue-white-red)
  if (is.null(palette)) {
    colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
                "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
    if (reverse) colors <- rev(colors)
    return(colors)
  }

  # paletteer palette
  if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("To use paletteer, install it with: install.packages('paletteer')")
    }

    # Try it as a discrete palette first
    colors <- tryCatch({
      raw_pal <- as.character(paletteer::paletteer_d(palette, n = n))
      vapply(raw_pal, .normalize_hex, character(1), USE.NAMES = FALSE)
    }, error = function(e) {
      # Then try it as a continuous palette
      tryCatch({
        raw_pal <- as.character(paletteer::paletteer_c(palette, n = n))
        vapply(raw_pal, .normalize_hex, character(1), USE.NAMES = FALSE)
      }, error = function(e2) {
        stop("Could not load palette '", palette, "': ", e2$message)
      })
    })

    if (reverse) colors <- rev(colors)
    return(colors)
  }

  # RColorBrewer palette
  if (is.character(palette) && length(palette) == 1 && startsWith(palette, "brewer:")) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("To use RColorBrewer, install it with: install.packages('RColorBrewer')")
    }

    nm <- sub("^brewer:", "", palette)
    colors <- tryCatch({
      maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
      RColorBrewer::brewer.pal(min(n, maxc), nm)
    }, error = function(e) {
      stop("Could not load brewer palette '", nm, "': ", e$message)
    })

    if (reverse) colors <- rev(colors)
    return(colors)
  }

  # Plain vector of colours
  if (is.character(palette) && length(palette) > 1) {
    if (reverse) palette <- rev(palette)
    return(palette)
  }

  # Bare palette name (try RColorBrewer)
  if (is.character(palette) && length(palette) == 1) {
    if (requireNamespace("RColorBrewer", quietly = TRUE) &&
        palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      maxc <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      colors <- RColorBrewer::brewer.pal(min(n, maxc), palette)
      if (reverse) colors <- rev(colors)
      return(colors)
    }
  }

  # Fallback default
  colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
              "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
  if (reverse) colors <- rev(colors)
  colors
}


#' Build the Colour Palette for Categorical Annotations
#'
#' @param levels Vector of levels (categories)
#' @param palette Palette specification (same format as get_heatmap_palette)
#'
#' @return Named vector of colours
get_annotation_palette <- function(levels, palette = NULL) {

  n <- length(levels)

  # Default palette
  default_palette <- c(
    "#457B9D", "#E63946", "#2A9D8F", "#E9C46A",
    "#9B5DE5", "#F4A261", "#264653", "#00BBF9",
    "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"
  )

  if (is.null(palette)) {
    pal <- grDevices::hcl.colors(n, "Dark 3")
    return(stats::setNames(pal, levels))
  }

  # paletteer palette
  if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("To use paletteer, install it with: install.packages('paletteer')")
    }
    colors <- tryCatch({
      raw_pal <- as.character(paletteer::paletteer_d(palette))
      vapply(raw_pal, .normalize_hex, character(1), USE.NAMES = FALSE)
    }, error = function(e) {
      stop("Could not load palette '", palette, "': ", e$message)
    })
    if (length(colors) < n) {
      colors <- rep(colors, length.out = n)
    }
    return(stats::setNames(colors[seq_len(n)], levels))
  }

  # RColorBrewer palette
  if (is.character(palette) && length(palette) == 1 && startsWith(palette, "brewer:")) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("To use RColorBrewer, install it with: install.packages('RColorBrewer')")
    }
    nm <- sub("^brewer:", "", palette)
    colors <- tryCatch({
      maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
      RColorBrewer::brewer.pal(max(3, min(n, maxc)), nm)
    }, error = function(e) {
      stop("Could not load brewer palette '", nm, "': ", e$message)
    })
    if (length(colors) < n) {
      colors <- rep(colors, length.out = n)
    }
    return(stats::setNames(colors[seq_len(n)], levels))
  }

  # Named vector of colours
  if (is.character(palette) && !is.null(names(palette))) {
    miss <- setdiff(levels, names(palette))
    if (length(miss) > 0) {
      stop("Missing colours for levels: ", paste(miss, collapse = ", "))
    }
    return(palette[levels])
  }

  # Unnamed vector of colours
  if (is.character(palette) && length(palette) >= 1) {
    if (length(palette) < n) {
      palette <- rep(palette, length.out = n)
    }
    return(stats::setNames(palette[seq_len(n)], levels))
  }

  # Fallback default
  if (length(default_palette) < n) {
    default_palette <- rep(default_palette, length.out = n)
  }
  stats::setNames(default_palette[seq_len(n)], levels)
}


# -----------------------------------------------------------------------------
# Main function: prepare the data for the heatmap
# -----------------------------------------------------------------------------

#' Prepare Long-Format Data for tidyHeatmap
#'
#' @param data Data frame in long format with columns:
#'   - SampleID: sample identifier
#'   - FeatureID: protein/feature identifier
#'   - Intensity: intensity value (log2)
#'   - Condition: experimental condition
#'   - Replicate: replicate number
#'   - sig_any: logical flagging significance (for mode="any")
#'   - adjP_*: adjusted p-value columns (for mode="target")
#' @param mode Filtering mode: "all", "any", or "target"
#' @param alpha Significance threshold for mode "target" (default: 0.05)
#' @param comparison Name of the comparison for mode "target" (e.g. "B-A")
#' @param feature_ids Vector of specific FeatureIDs (when supplied, mode/alpha are ignored)
#' @param scale_data Scaling type: "none", "row", "column" (default: "row")
#' @param sample_order Sample order: "clustering", "condition", or a custom vector
#' @param condition_order Condition order when sample_order = "condition"
#'
#' @return Data frame in long format, ready for tidyHeatmap
prepare_heatmap_data <- function(data,
                                 mode = c("all", "any", "target"),
                                 alpha = 0.05,
                                 comparison = NULL,
                                 feature_ids = NULL,
                                 scale_data = c("row", "none", "column"),
                                 sample_order = c("clustering", "condition", "custom"),
                                 condition_order = NULL) {

  mode <- match.arg(mode)
  scale_data <- match.arg(scale_data)

  if (is.character(sample_order) && length(sample_order) == 1) {
    sample_order <- match.arg(sample_order, c("clustering", "condition", "custom"))
  }

  # Validate the required columns
  required_cols <- c("SampleID", "FeatureID", "Intensity", "Condition", "Replicate")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  data <- as.data.frame(data)

  # For mode = "target", keep only the samples belonging to the compared conditions
  if (mode == "target" && !is.null(comparison)) {
    # Extract the conditions from the comparison (e.g. "B-A" -> c("B", "A"))
    conds <- unique(trimws(strsplit(comparison, "[-|:]")[[1]]))
    if (length(conds) >= 2) {
      data <- data[data$Condition %in% conds, , drop = FALSE]
      # Update condition_order so it only lists the relevant conditions
      if (!is.null(condition_order)) {
        condition_order <- condition_order[condition_order %in% conds]
      } else {
        condition_order <- conds
      }
    }
  }

  # Get the feature IDs according to the mode, or use the ones supplied
  if (!is.null(feature_ids) && length(feature_ids) > 0) {
    # Use the supplied IDs directly
    available_ids <- unique(data$FeatureID)
    ids <- feature_ids[feature_ids %in% available_ids]
    if (length(ids) < length(feature_ids)) {
      missing <- setdiff(feature_ids, available_ids)
      warning("FeatureIDs not found in the data (ignored): ",
              paste(head(missing, 5), collapse = ", "),
              if (length(missing) > 5) paste0(" ... and ", length(missing) - 5, " more"))
    }
  } else {
    # Use the mode-based filtering
    ids <- .get_feature_ids(data, mode = mode, alpha = alpha, comparison = comparison)
  }

  if (length(ids) < 2) {
    stop("Subset does not have enough proteins (minimum 2). Found: ", length(ids))
  }

  # Filter the data
  dt <- data[data$FeatureID %in% ids, , drop = FALSE]
  dt <- dt[is.finite(dt$Intensity) & !is.na(dt$Intensity), , drop = FALSE]


  # Build the long-format data frame for the heatmap
  # Include adjP when mode = "target", for the row annotation
  if (mode == "target" && !is.null(comparison)) {
    adjp_colname <- .adjp_col(comparison)
    if (adjp_colname %in% names(dt)) {
      hm_data <- dt %>%
        select(SampleID, FeatureID, Intensity, Condition, Replicate, all_of(adjp_colname)) %>%
        rename(adjP = all_of(adjp_colname)) %>%
        distinct()
    } else {
      hm_data <- dt %>%
        select(SampleID, FeatureID, Intensity, Condition, Replicate) %>%
        distinct()
    }
  } else {
    hm_data <- dt %>%
      select(SampleID, FeatureID, Intensity, Condition, Replicate) %>%
      distinct()
  }

  # Order the samples according to the requested criterion
  if (is.character(sample_order) && length(sample_order) == 1 && sample_order == "condition") {
    # Sort by Condition and then by Replicate
    if (!is.null(condition_order)) {
      hm_data <- hm_data %>%
        mutate(Condition = factor(Condition, levels = condition_order))
    }
    hm_data <- hm_data %>%
      arrange(Condition, Replicate) %>%
      mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
  } else if (is.character(sample_order) && length(sample_order) > 1) {
    # Custom order (vector of SampleIDs)
    # Keep only the specified samples
    samples_in_data <- unique(hm_data$SampleID)
    valid_samples <- sample_order[sample_order %in% samples_in_data]

    if (length(valid_samples) == 0) {
      stop("None of the samples given in sample_order is present in the data")
    }

    if (!all(sample_order %in% samples_in_data)) {
      missing <- setdiff(sample_order, samples_in_data)
      warning("Samples not found in the data (ignored): ", paste(missing, collapse = ", "))
    }

    # Subset the data down to the specified samples
    hm_data <- hm_data %>%
      filter(SampleID %in% valid_samples) %>%
      mutate(SampleID = factor(SampleID, levels = valid_samples))
  }
  # When sample_order == "clustering" we let tidyHeatmap do the clustering

  # Apply the scaling if requested
  if (scale_data != "none") {
    # Pivot to a matrix in order to scale (use mean to collapse possible duplicates)
    mat_wide <- hm_data %>%
      select(SampleID, FeatureID, Intensity) %>%
      pivot_wider(names_from = SampleID, values_from = Intensity, values_fn = mean) %>%
      as.data.frame()

    rownames_feat <- mat_wide$FeatureID
    mat <- as.matrix(mat_wide[, -1])
    rownames(mat) <- rownames_feat

    # NA-aware scaling: base scale() does NOT accept na.rm, so a single NA cell
    # would turn the whole row/column into NA. Centre/scale by hand instead.
    if (scale_data == "row") {
      ctr <- rowMeans(mat, na.rm = TRUE)
      sdv <- apply(mat, 1, sd, na.rm = TRUE)
      # Constant features (sd 0/NA) -> undefined z-score, and they would break
      # the row clustering (all-NA row). Drop them before scaling.
      keep <- is.finite(sdv) & sdv > 0
      if (any(!keep)) {
        message(sprintf(
          "Heatmap: %d features with 0/NA variance dropped before the z-score.",
          sum(!keep)))
        mat <- mat[keep, , drop = FALSE]; ctr <- ctr[keep]; sdv <- sdv[keep]
      }
      mat <- sweep(mat, 1, ctr, "-")
      mat <- sweep(mat, 1, sdv, "/")
    } else if (scale_data == "column") {
      ctr <- colMeans(mat, na.rm = TRUE)
      sdv <- apply(mat, 2, sd, na.rm = TRUE)
      sdv[!is.finite(sdv) | sdv == 0] <- NA
      mat <- sweep(mat, 2, ctr, "-")
      mat <- sweep(mat, 2, sdv, "/")
    }

    # Replace NaN with NA
    mat[is.nan(mat)] <- NA

    # Back to long format
    mat_df <- as.data.frame(mat)
    mat_df$FeatureID <- rownames(mat)

    scaled_long <- mat_df %>%
      pivot_longer(cols = -FeatureID, names_to = "SampleID", values_to = "Intensity")

    # Merge the metadata back in (including adjP when present)
    if ("adjP" %in% names(hm_data)) {
      metadata <- hm_data %>%
        select(SampleID, FeatureID, Condition, Replicate, adjP) %>%
        distinct()

      hm_data <- scaled_long %>%
        left_join(metadata, by = c("SampleID", "FeatureID"))
    } else {
      metadata <- hm_data %>%
        select(SampleID, Condition, Replicate) %>%
        distinct()

      hm_data <- scaled_long %>%
        left_join(metadata, by = "SampleID")
    }

    # Restore the factor ordering where applicable
    if (is.character(sample_order) && length(sample_order) == 1 && sample_order == "condition") {
      if (!is.null(condition_order)) {
        hm_data <- hm_data %>%
          mutate(Condition = factor(Condition, levels = condition_order))
      }
      hm_data <- hm_data %>%
        arrange(Condition, as.numeric(Replicate)) %>%
        mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
    } else if (is.character(sample_order) && length(sample_order) > 1) {
      # Restore the custom order
      hm_data <- hm_data %>%
        mutate(SampleID = factor(SampleID, levels = sample_order))
    }
  }

  hm_data
}


# -----------------------------------------------------------------------------
# Main function: heatmap built with tidyHeatmap
# -----------------------------------------------------------------------------

#' Create a tidyHeatmap Heatmap for Proteomics
#'
#' @param data Data frame in long format (see prepare_heatmap_data for the structure)
#' @param mode Protein filtering mode: "all", "any", or "target"
#' @param alpha Significance threshold for DEP proteins (default: 0.05)
#' @param comparison Name of the comparison for mode "target" (e.g. "B-A")
#' @param feature_ids Vector of specific FeatureIDs to show (default: NULL, uses the mode-based filtering).
#'   When supplied, only these proteins are shown and the mode/alpha filtering is ignored.
#' @param row_annotation Row annotations (FeatureIDs). It can be:
#'   - Path to a TSV file with a "FeatureID" column plus additional categorical columns
#'   - Data frame with the same structure
#'   - NULL: no row annotations (default)
#' @param row_annotation_cols Vector of column names to display as annotations.
#'   Default: NULL (uses every column except FeatureID)
#' @param row_annotation_palette Named list of palettes, one per annotation.
#'   Example: list(Pathway = "brewer:Set1", Function = c("red", "blue", "green"))
#' @param row_annotation_size Width of the row annotation bars. It can be:
#'   - Number: interpreted as centimetres (e.g. 0.3 = 0.3cm)
#'   - A unit object: grid::unit(0.3, "cm")
#'   - NULL: use tidyHeatmap's default value
#' @param row_annotation_name_size Font size of the row annotation names (default: 8)
#' @param row_order_by Row order. It can be:
#'   - NULL: no particular ordering (default)
#'   - "clustering": order by hierarchical clustering
#'   - A column name: order by that annotation (e.g. "Specie")
#'   - A vector of columns: order sequentially (e.g. c("Specie", "Process", "Function"))
#'   - Include "adjP" when mode = "target" to order by adjusted p-value
#' @param split_rows_by Name of the annotation column used to split the rows into groups (default: NULL)
#' @param scale_data Scaling type: "none", "row", "column" (default: "row")
#' @param sample_order Sample order: "clustering", "condition", or a custom vector
#' @param condition_order Condition order when sample_order = "condition"
#' @param cluster_rows Cluster the rows (default: FALSE)
#' @param cluster_columns Cluster the columns (default: TRUE)
#' @param show_row_names Show the row names (default: TRUE when <= 50 proteins)
#' @param show_column_names Show the column names (default: TRUE)
#' @param palette_value Palette for the heatmap values (see get_heatmap_palette)
#' @param palette_annotation Palette for the Condition annotation
#' @param reverse_palette Reverse the value palette (default: FALSE)
#' @param row_title Title for the rows
#' @param column_title Title for the columns
#' @param row_title_size Font size of the row title (default: 10)
#' @param column_title_size Font size of the column title (default: 10)
#' @param show_row_title Show the row title (default: TRUE)
#' @param show_column_title Show the column title (default: TRUE)
#' @param show_annotation Show the Condition annotation (default: TRUE)
#' @param split_by_condition Split the heatmap by condition (default: FALSE)
#' @param show_adjp_annotation Show the adjP row annotation when mode="target" (default: TRUE)
#' @param palette_adjp Palette for the adjP annotation. It can be:
#'   - NULL: use the default palette (red-orange-white)
#'   - A vector of 3 colours: c(color_0, color_middle, color_0.05)
#'   - String "brewer:name": use an RColorBrewer palette
#'   - String "package::palette": use a paletteer palette
#' @param row_names_size Font size of the row names (default: 7)
#' @param column_names_size Font size of the column names (default: 9)
#' @param column_names_rotation Rotation of the column names in degrees (default: 45)
#' @param column_dend_height Height of the column dendrogram. It can be:
#'   - Number: interpreted as millimetres (e.g. 30 = 30mm)
#'   - A unit object: grid::unit(2, "cm")
#'   - NULL: use ComplexHeatmap's default value
#' @param row_dend_width Width of the row dendrogram. Same format as column_dend_height
#' @param show_heatmap_legend Show the heatmap legend (default: TRUE)
#' @param show_annotation_legend Show the annotation legends (default: TRUE)
#' @param border_color Border colour of the heatmap cells. It can be:
#'   - NULL or FALSE: no border (default)
#'   - TRUE: black border
#'   - A colour string: that specific colour (e.g. "black", "grey", "#CCCCCC")
#' @param heatmap_title Main heatmap title (default: NULL, no title)
#' @param heatmap_title_size Font size of the main title (default: 14)
#' @param heatmap_title_face Font face of the title: "plain", "bold", "italic", "bold.italic" (default: "bold")
#' @param export_path Path used to export the heatmap data to TSV (default: NULL, nothing is exported).
#'   The file contains FeatureID, the per-sample intensity values, and the metadata (adjP where applicable).
#' @param export_file Path used to export the plot. The format is taken from the extension:
#'   - .png: PNG image (raster)
#'   - .svg: SVG image (vector)
#'   - .pdf: PDF document (vector)
#' @param plot_width Plot width in inches (default: 10).
#'   To get pixels: pixels = inches x dpi (e.g. 10" x 300dpi = 3000px)
#' @param plot_height Plot height in inches (default: 8).
#'   To get pixels: pixels = inches x dpi (e.g. 8" x 300dpi = 2400px)
#' @param export_dpi Resolution for PNG output in dots per inch (default: 300).
#'   Higher dpi = more detail. Common values: 72 (web), 150 (draft), 300 (publication)
#'
#' @return tidyHeatmap/ComplexHeatmap object
#'
#' @examples
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("tidyHeatmap", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#'   # Proteins significant in the B-A comparison, samples grouped by condition
#'   hm <- proteomics_heatmap(
#'     res$PCA_Input,
#'     mode            = "target",
#'     comparison      = "B-A",
#'     scale_data      = "row",
#'     sample_order    = "condition",
#'     condition_order = c("A", "B", "D")
#'   )
#'   print(class(hm))   # print(hm) itself draws it on the current device
#' }
#'
#' @export
proteomics_heatmap <- function(data,
                               mode = c("all", "any", "target"),
                               alpha = 0.05,
                               comparison = NULL,
                               feature_ids = NULL,
                               row_annotation = NULL,
                               row_annotation_cols = NULL,
                               row_annotation_palette = NULL,
                               row_annotation_size = NULL,
                               row_annotation_name_size = 8,
                               row_order_by = NULL,
                               split_rows_by = NULL,
                               scale_data = c("row", "none", "column"),
                               sample_order = "clustering",
                               condition_order = NULL,
                               cluster_rows = FALSE,
                               cluster_columns = TRUE,
                               show_row_names = NULL,
                               show_column_names = TRUE,
                               palette_value = NULL,
                               palette_annotation = NULL,
                               reverse_palette = FALSE,
                               row_title = "Proteins",
                               column_title = "Samples",
                               row_title_size = 10,
                               column_title_size = 10,
                               show_row_title = TRUE,
                               show_column_title = TRUE,
                               show_annotation = TRUE,
                               split_by_condition = FALSE,
                               show_adjp_annotation = TRUE,
                               palette_adjp = NULL,
                               row_names_size = 7,
                               column_names_size = 9,
                               column_names_rotation = 45,
                               column_dend_height = NULL,
                               row_dend_width = NULL,
                               show_heatmap_legend = TRUE,
                               show_annotation_legend = TRUE,
                               border_color = NULL,
                               heatmap_title = NULL,
                               heatmap_title_size = 14,
                               heatmap_title_face = c("bold", "plain", "italic", "bold.italic"),
                               export_path = NULL,
                               export_file = NULL,
                               plot_width = 10,
                               plot_height = 8,
                               export_dpi = 300) {

  mode <- match.arg(mode)
  scale_data <- match.arg(scale_data)
  heatmap_title_face <- match.arg(heatmap_title_face)

  # Silence ComplexHeatmap's messages (use_raster, magick)
  # Set for the whole session, since the messages appear when drawing, not when creating
  ComplexHeatmap::ht_opt(message = FALSE)

  # Resolve border_color
  rect_gp <- NULL
  if (!is.null(border_color) && !isFALSE(border_color)) {
    if (isTRUE(border_color)) {
      rect_gp <- grid::gpar(col = "black")
    } else {
      rect_gp <- grid::gpar(col = border_color)
    }
  }

  # Resolve the titles (hide them when show_*_title is FALSE)
  if (!show_row_title) {
    row_title <- NULL
  }
  if (!show_column_title) {
    column_title <- NULL
  }

  # Convert the dendrogram sizes to grid units when they are given as numbers
  if (!is.null(column_dend_height) && is.numeric(column_dend_height)) {
    column_dend_height <- grid::unit(column_dend_height, "mm")
  }
  if (!is.null(row_dend_width) && is.numeric(row_dend_width)) {
    row_dend_width <- grid::unit(row_dend_width, "mm")
  }

  # ---------------------------------------------------------------------------
  # 1) Prepare the data
  # ---------------------------------------------------------------------------

  hm_data <- prepare_heatmap_data(
    data = data,
    mode = mode,
    alpha = alpha,
    comparison = comparison,
    feature_ids = feature_ids,
    scale_data = scale_data,
    sample_order = sample_order,
    condition_order = condition_order
  )

  # ---------------------------------------------------------------------------
  # 1b) Export the data to TSV if requested
  # ---------------------------------------------------------------------------

  if (!is.null(export_path) && nzchar(export_path)) {
    # Build a wide matrix with FeatureID as rows and SampleID as columns
    export_wide <- hm_data %>%
      select(FeatureID, SampleID, Intensity) %>%
      tidyr::pivot_wider(
        names_from = SampleID,
        values_from = Intensity,
        values_fn = mean
      )

    # Add adjP when present
    if ("adjP" %in% names(hm_data)) {
      adjp_data <- hm_data %>%
        select(FeatureID, adjP) %>%
        distinct()
      export_wide <- export_wide %>%
        left_join(adjp_data, by = "FeatureID")
    }

    # Write the TSV
    readr::write_tsv(export_wide, export_path)
    message("Data exported to: ", export_path)
  }

  # ---------------------------------------------------------------------------
  # 1c) Process the row annotations
  # ---------------------------------------------------------------------------

  row_annot_data <- NULL
  row_annot_colors <- list()
  row_split_vector <- NULL

  if (!is.null(row_annotation)) {
    # Load the annotations when a file path is given
    if (is.character(row_annotation) && length(row_annotation) == 1 && file.exists(row_annotation)) {
      row_annot_data <- readr::read_tsv(row_annotation, show_col_types = FALSE)
    } else if (is.data.frame(row_annotation)) {
      row_annot_data <- as.data.frame(row_annotation)
    } else {
      warning("row_annotation must be a path to a TSV file or a data.frame")
    }

    if (!is.null(row_annot_data)) {
      # Validate the FeatureID column
      if (!("FeatureID" %in% names(row_annot_data))) {
        stop("The annotation file must have a 'FeatureID' column")
      }

      # Work out which columns to use
      all_annot_cols <- setdiff(names(row_annot_data), "FeatureID")
      if (is.null(row_annotation_cols)) {
        row_annotation_cols <- all_annot_cols
      } else {
        missing_cols <- setdiff(row_annotation_cols, all_annot_cols)
        if (length(missing_cols) > 0) {
          warning("Annotation columns not found: ", paste(missing_cols, collapse = ", "))
          row_annotation_cols <- intersect(row_annotation_cols, all_annot_cols)
        }
      }

      # Turn empty values into the "NA" string (not NA_character_) so they can be coloured white
      for (col in row_annotation_cols) {
        values <- row_annot_data[[col]]
        if (is.character(values)) {
          # Convert empty strings, whitespace and real NAs into the "NA" string
          row_annot_data[[col]] <- ifelse(
            is.na(values) | trimws(values) == "",
            "NA",
            values
          )
        }
      }

      # Keep only the FeatureIDs present in the data
      feature_ids_in_data <- unique(hm_data$FeatureID)
      row_annot_data <- row_annot_data[row_annot_data$FeatureID %in% feature_ids_in_data, , drop = FALSE]

      # Build a colour palette for each annotation (mapping "NA" to white)
      for (col in row_annotation_cols) {
        # Get the unique levels
        all_levels <- unique(row_annot_data[[col]])
        has_na <- "NA" %in% all_levels

        # Separate the real levels from "NA" and order them: real ones first, "NA" last
        real_levels <- sort(all_levels[all_levels != "NA"])
        if (has_na) {
          ordered_levels <- c(real_levels, "NA")
        } else {
          ordered_levels <- real_levels
        }

        # Turn the column into a factor with the ordered levels
        row_annot_data[[col]] <- factor(row_annot_data[[col]], levels = ordered_levels)

        # Get the colours for the real levels
        if (!is.null(row_annotation_palette) && col %in% names(row_annotation_palette)) {
          real_colors <- get_annotation_palette(real_levels, row_annotation_palette[[col]])
        } else {
          real_colors <- get_annotation_palette(real_levels, NULL)
        }

        # Build the final palette in the EXACT order of the factor levels;
        # this is what makes tidyHeatmap assign the colours correctly
        final_palette <- character(length(ordered_levels))
        names(final_palette) <- ordered_levels
        for (lvl in ordered_levels) {
          if (lvl == "NA") {
            final_palette[lvl] <- "#FFFFFF"
          } else {
            final_palette[lvl] <- real_colors[lvl]
          }
        }
        row_annot_colors[[col]] <- final_palette
      }

      # Join the annotations onto hm_data
      cols_to_join <- c("FeatureID", row_annotation_cols)
      hm_data <- hm_data %>%
        left_join(row_annot_data[, cols_to_join, drop = FALSE], by = "FeatureID")

      # Make sure the factors survive the join
      for (col in row_annotation_cols) {
        if (col %in% names(hm_data) && col %in% names(row_annot_data)) {
          hm_data[[col]] <- factor(hm_data[[col]], levels = levels(row_annot_data[[col]]))
        }
      }

      # Order the rows by annotation(s) when requested (except "clustering")
      if (!is.null(row_order_by) && !identical(row_order_by, "clustering")) {
        # Work out which columns are valid to order by:
        # the annotation columns plus adjP (when present in hm_data)
        valid_order_cols <- row_annotation_cols
        if ("adjP" %in% names(hm_data)) {
          valid_order_cols <- c(valid_order_cols, "adjP")
        }

        # Keep only the valid columns from the row_order_by vector
        order_cols <- row_order_by[row_order_by %in% valid_order_cols]

        if (length(order_cols) > 0) {
          # Assemble the data to sort on (annotations plus adjP where needed)
          order_data <- row_annot_data

          # Add adjP to order_data when it is needed for the ordering
          if ("adjP" %in% order_cols && "adjP" %in% names(hm_data)) {
            adjp_data <- hm_data %>%
              select(FeatureID, adjP) %>%
              distinct()
            order_data <- order_data %>%
              left_join(adjp_data, by = "FeatureID")
          }

          # Sort by several columns sequentially
          order_syms <- rlang::syms(order_cols)
          order_data_sorted <- order_data %>%
            arrange(!!!order_syms)

          # Get the FeatureID order
          feature_order <- unique(order_data_sorted$FeatureID)

          # Turn FeatureID into a factor carrying the right order
          hm_data <- hm_data %>%
            mutate(FeatureID = factor(FeatureID, levels = feature_order)) %>%
            arrange(FeatureID)  # Sort the data explicitly
        }
      }

      # Build row_split when requested (AFTER the ordering)
      if (!is.null(split_rows_by) && split_rows_by %in% row_annotation_cols) {
        # Get the unique FeatureIDs in the current order of hm_data
        if (!is.null(row_order_by)) {
          # When an order was applied, use the factor levels
          ordered_features <- levels(hm_data$FeatureID)
        } else {
          ordered_features <- unique(as.character(hm_data$FeatureID))
        }

        # Build the FeatureID -> split value mapping
        feature_to_split <- row_annot_data %>%
          select(FeatureID, all_of(split_rows_by)) %>%
          distinct()
        feature_to_split <- stats::setNames(
          feature_to_split[[split_rows_by]],
          feature_to_split$FeatureID
        )

        # Build the split vector in the right order
        split_values <- feature_to_split[ordered_features]

        # Determine the order of the split levels (as they appear in the sorted data)
        split_levels <- unique(split_values)
        split_levels <- split_levels[!is.na(split_levels)]

        row_split_vector <- factor(split_values, levels = split_levels)
      }
    }
  }

  # Order by adjP when requested and there are no row annotations
  # (when there are annotations, the ordering is already handled above)
  if (!is.null(row_order_by) && !identical(row_order_by, "clustering") &&
      (is.null(row_annot_data) || length(row_annotation_cols) == 0)) {
    # Check whether ordering by adjP was requested
    if ("adjP" %in% row_order_by && "adjP" %in% names(hm_data)) {
      # Get the FeatureID order by adjP
      adjp_order <- hm_data %>%
        select(FeatureID, adjP) %>%
        distinct() %>%
        arrange(adjP)
      feature_order <- adjp_order$FeatureID

      # Turn FeatureID into a factor carrying the right order
      hm_data <- hm_data %>%
        mutate(FeatureID = factor(FeatureID, levels = feature_order)) %>%
        arrange(FeatureID)
    }
  }

  # Count the proteins, used for the auto-configuration below
  n_proteins <- length(unique(hm_data$FeatureID))
  n_samples <- length(unique(hm_data$SampleID))

  # Decide automatically whether to show the row names
  if (is.null(show_row_names)) {
    show_row_names <- n_proteins <= 50
  }

  # ---------------------------------------------------------------------------
  # 2) Build the colour palettes
  # ---------------------------------------------------------------------------

  # Palette for the values
  value_colors <- get_heatmap_palette(
    palette = palette_value,
    n = 11,
    reverse = reverse_palette
  )

  # Build a colorRamp2 function for the values
  if (requireNamespace("circlize", quietly = TRUE)) {
    # Determine the value range
    val_range <- range(hm_data$Intensity, na.rm = TRUE)

    if (scale_data != "none") {
      # For scaled data, use a symmetric range
      abs_max <- max(abs(val_range), na.rm = TRUE)
      val_breaks <- seq(-abs_max, abs_max, length.out = length(value_colors))
    } else {
      val_breaks <- seq(val_range[1], val_range[2], length.out = length(value_colors))
    }

    palette_func <- circlize::colorRamp2(val_breaks, value_colors)
  } else {
    palette_func <- value_colors
  }

  # Palette for the Condition annotation
  condition_levels <- unique(hm_data$Condition)
  if (!is.null(condition_order)) {
    condition_levels <- condition_order[condition_order %in% condition_levels]
  }
  annotation_colors <- get_annotation_palette(condition_levels, palette_annotation)

  # ---------------------------------------------------------------------------
  # 3) Configure the clustering
  # ---------------------------------------------------------------------------

  # Decide whether to cluster the columns
  cluster_cols_final <- cluster_columns
  if (is.character(sample_order) && (length(sample_order) > 1 || sample_order != "clustering")) {
    cluster_cols_final <- FALSE
  }

  # Decide whether to cluster the rows
  cluster_rows_final <- cluster_rows
  if (!is.null(row_order_by)) {
    if (identical(row_order_by, "clustering")) {
      # row_order_by = "clustering" turns on row clustering
      cluster_rows_final <- TRUE
    } else {
      # Ordering by specific columns turns clustering off
      cluster_rows_final <- FALSE
    }
  }
  # split_rows_by is compatible with clustering (it clusters within each split)

  # ---------------------------------------------------------------------------
  # 4) Build the heatmap with tidyHeatmap
  # ---------------------------------------------------------------------------

  # Convert to a tibble holding base R types (avoids trouble with arrow/parquet classes):
  # tidyHeatmap uses `class(x) %in% ...`, which fails when class() returns several values
  hm_data <- as.data.frame(lapply(hm_data, function(col) {
    if (is.factor(col)) return(factor(as.character(col), levels = levels(col)))
    if (is.logical(col)) return(as.logical(col))
    if (is.numeric(col)) return(as.numeric(col))
    if (is.character(col)) return(as.character(col))
    col
  }), stringsAsFactors = FALSE)
  hm_data <- tibble::as_tibble(hm_data)

  # Build column_split when requested (to visually separate the conditions)
  col_split_vector <- NULL
  if (split_by_condition) {
    # Get the unique samples in the order in which they appear in the data
    sample_info <- hm_data %>%
      select(SampleID, Condition) %>%
      distinct()

    # Preserve the original SampleID order of the data
    sample_order_vec <- unique(as.character(hm_data$SampleID))
    sample_info <- sample_info[match(sample_order_vec, as.character(sample_info$SampleID)), ]
    col_split_vector <- factor(sample_info$Condition, levels = unique(sample_info$Condition))
  }

  # Assemble the optional arguments
  extra_args <- list()
  if (!is.null(column_dend_height)) {
    extra_args$column_dend_height <- column_dend_height
  }
  if (!is.null(row_dend_width)) {
    extra_args$row_dend_width <- row_dend_width
  }
  if (!is.null(rect_gp)) {
    extra_args$rect_gp <- rect_gp
  }
  if (!is.null(row_split_vector)) {
    extra_args$row_split <- row_split_vector
  }

  # Build the base heatmap
  hm_args <- c(
    list(
      .data = hm_data,
      .row = rlang::sym("FeatureID"),
      .column = rlang::sym("SampleID"),
      .value = rlang::sym("Intensity"),
      scale = "none",
      cluster_rows = cluster_rows_final,
      cluster_columns = cluster_cols_final,
      palette_value = palette_func,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_gp = grid::gpar(fontsize = row_names_size),
      column_names_gp = grid::gpar(fontsize = column_names_size),
      column_names_rot = column_names_rotation,
      row_title = row_title,
      column_title = column_title,
      row_title_gp = grid::gpar(fontsize = row_title_size),
      column_title_gp = grid::gpar(fontsize = column_title_size),
      column_split = col_split_vector,
      show_heatmap_legend = show_heatmap_legend
    ),
    extra_args
  )

  hm <- do.call(tidyHeatmap::heatmap, hm_args)

  # Add the Condition annotation as a colour bar when requested
  if (show_annotation) {
    hm <- hm %>%
      tidyHeatmap::annotation_tile(
        Condition,
        palette = annotation_colors,
        show_legend = show_annotation_legend
      )
  }

  # Add the custom row annotations (suppressWarnings so NAs are handled quietly)
  if (!is.null(row_annot_data) && length(row_annotation_cols) > 0) {
    # Turn size into a unit when given as a number
    annot_size <- NULL
    if (!is.null(row_annotation_size)) {
      if (is.numeric(row_annotation_size)) {
        annot_size <- grid::unit(row_annotation_size, "cm")
      } else {
        annot_size <- row_annotation_size
      }
    }

    # Assemble annotation_name_gp
    annot_name_gp <- grid::gpar(fontsize = row_annotation_name_size)

    # Identify the valid columns
    valid_cols <- row_annotation_cols[row_annotation_cols %in% names(hm_data)]
    n_cols <- length(valid_cols)

    for (i in seq_along(valid_cols)) {
      col <- valid_cols[i]
      # Only set size on the last annotation (avoids a tidyHeatmap warning)
      use_size <- if (i == n_cols) annot_size else NULL

      hm <- hm %>%
        tidyHeatmap::annotation_tile(
          !!rlang::sym(col),
          palette = row_annot_colors[[col]],
          size = use_size,
          annotation_name_gp = annot_name_gp,
          show_legend = show_annotation_legend
        )
    }
  }

  # Add the adjP annotation for mode = "target" (a row annotation)
  if (show_adjp_annotation && mode == "target" && "adjP" %in% names(hm_data)) {
    # Get the colours for the adjP palette
    adjp_colors <- c("#67001F", "#F4A582", "#F7F7F7")  # Default: dark red -> orange -> white

    if (!is.null(palette_adjp)) {
      if (is.character(palette_adjp) && length(palette_adjp) >= 3) {
        # Custom vector of colours
        adjp_colors <- palette_adjp[seq_len(3)]
      } else if (is.character(palette_adjp) && length(palette_adjp) == 1) {
        # paletteer or RColorBrewer palette
        if (grepl("::", palette_adjp)) {
          # Paletteer
          if (requireNamespace("paletteer", quietly = TRUE)) {
            raw_pal <- tryCatch({
              as.character(paletteer::paletteer_c(palette_adjp, n = 3))
            }, error = function(e) {
              tryCatch({
                as.character(paletteer::paletteer_d(palette_adjp, n = 3))
              }, error = function(e2) NULL)
            })
            if (!is.null(raw_pal) && length(raw_pal) >= 3) {
              adjp_colors <- vapply(raw_pal[seq_len(3)], .normalize_hex,
                                    character(1), USE.NAMES = FALSE)
            }
          }
        } else if (startsWith(palette_adjp, "brewer:")) {
          # RColorBrewer
          if (requireNamespace("RColorBrewer", quietly = TRUE)) {
            nm <- sub("^brewer:", "", palette_adjp)
            raw_pal <- tryCatch({
              RColorBrewer::brewer.pal(3, nm)
            }, error = function(e) NULL)
            if (!is.null(raw_pal)) {
              adjp_colors <- raw_pal
            }
          }
        }
      }
    }

    # Build the palette for the p-values (low values = more significant)
    adjp_palette <- circlize::colorRamp2(
      c(0, 0.01, 0.05),
      adjp_colors
    )

    hm <- hm %>%
      tidyHeatmap::annotation_tile(
        adjP,
        palette = adjp_palette,
        show_legend = show_annotation_legend
      )
  }

  # ---------------------------------------------------------------------------
  # 5) Add the main title when one was given
  # ---------------------------------------------------------------------------

  if (!is.null(heatmap_title) && nzchar(heatmap_title)) {
    hm <- wrap_heatmap_with_title(
      hm = hm,
      title = heatmap_title,
      title_size = heatmap_title_size,
      title_face = heatmap_title_face
    )
  }

  # ---------------------------------------------------------------------------
  # 6) Export the plot when a path was given
  # ---------------------------------------------------------------------------

  if (!is.null(export_file) && nzchar(export_file)) {
    # Take the format from the extension
    file_ext <- tolower(tools::file_ext(export_file))

    # Create the directory if it does not exist
    export_dir <- dirname(export_file)
    if (!dir.exists(export_dir) && export_dir != ".") {
      dir.create(export_dir, recursive = TRUE)
    }

    # Compute the pixel size for PNG output (inches x dpi)
    width_pixels <- round(plot_width * export_dpi)
    height_pixels <- round(plot_height * export_dpi)

    # Helper that draws the heatmap according to its type
    draw_heatmap <- function(hm_obj) {
      if (inherits(hm_obj, "proteomics_heatmap")) {
        # Object carrying a title - the print method already handles everything
        print(hm_obj)
      } else if (inherits(hm_obj, c("Heatmap", "HeatmapList"))) {
        # Native ComplexHeatmap
        ComplexHeatmap::draw(hm_obj)
      } else if (inherits(hm_obj, "InputHeatmap")) {
        # tidyHeatmap - convert to ComplexHeatmap and draw
        ht <- tidyHeatmap::as_ComplexHeatmap(hm_obj)
        ComplexHeatmap::draw(ht)
      } else {
        # Fallback for any other type
        print(hm_obj)
      }
    }

    tryCatch({
      if (file_ext == "png") {
        grDevices::png(
          filename = export_file,
          width = width_pixels,
          height = height_pixels,
          res = export_dpi
        )
        draw_heatmap(hm)
        grDevices::dev.off()
        message("Heatmap exported to PNG: ", export_file,
                " (", width_pixels, "x", height_pixels, "px)")

      } else if (file_ext == "svg") {
        grDevices::svg(
          filename = export_file,
          width = plot_width,
          height = plot_height
        )
        draw_heatmap(hm)
        grDevices::dev.off()
        message("Heatmap exported to SVG: ", export_file,
                " (", plot_width, "x", plot_height, "in)")

      } else if (file_ext == "pdf") {
        grDevices::pdf(
          file = export_file,
          width = plot_width,
          height = plot_height
        )
        draw_heatmap(hm)
        grDevices::dev.off()
        message("Heatmap exported to PDF: ", export_file,
                " (", plot_width, "x", plot_height, "in)")

      } else {
        warning("Unsupported format: ", file_ext, ". Use .png, .svg or .pdf")
      }
    }, error = function(e) {
      # Make sure the graphics device is closed even if something fails
      try(grDevices::dev.off(), silent = TRUE)
      warning("Could not export the heatmap: ", e$message)
    })
  }

  hm
}


# -----------------------------------------------------------------------------
# Wrapper function: build a list of heatmaps
# -----------------------------------------------------------------------------

#' Build a List of Heatmaps for Several Subsets
#'
#' Automatically builds heatmaps for "all", "any" and/or specific comparisons.
#'
#' @param data Data frame in long format (see prepare_heatmap_data for the structure)
#' @param modes Vector of modes to generate: "all", "any", and/or comparison names
#'   (default: c("all", "any"))
#' @param alpha Significance threshold for DEP proteins (default: 0.05)
#' @param feature_ids Vector of specific FeatureIDs to show (default: NULL)
#' @param row_annotation Row annotations (TSV path or data.frame)
#' @param row_annotation_cols Columns to use as row annotations
#' @param row_annotation_palette List of palettes for the row annotations
#' @param row_annotation_size Width of the row annotation bars (a number in cm, or a unit)
#' @param row_annotation_name_size Font size of the row annotation names (default: 8)
#' @param row_order_by Row order: "clustering", annotation column(s), or include "adjP" for mode="target"
#' @param split_rows_by Annotation column used to split the rows into groups
#' @param scale_data Scaling type: "none", "row", "column" (default: "row")
#' @param sample_order Sample order: "clustering", "condition", or a custom vector
#' @param condition_order Condition order when sample_order = "condition"
#' @param cluster_rows Cluster the rows (default: FALSE)
#' @param cluster_columns Cluster the columns (default: TRUE)
#' @param show_row_names Show the row names (default: auto)
#' @param show_column_names Show the column names (default: TRUE)
#' @param palette_value Palette for the heatmap values
#' @param palette_annotation Palette for the Condition annotation
#' @param reverse_palette Reverse the value palette (default: FALSE)
#' @param row_title_size Font size of the row title (default: 10)
#' @param column_title_size Font size of the column title (default: 10)
#' @param show_row_title Show the row title (default: TRUE)
#' @param show_column_title Show the column title (default: TRUE)
#' @param show_annotation Show the Condition annotation (default: TRUE)
#' @param split_by_condition Split the heatmap by condition (default: FALSE)
#' @param show_adjp_annotation Show the adjP row annotation for comparisons (default: TRUE)
#' @param palette_adjp Palette for the adjP annotation (see proteomics_heatmap)
#' @param row_names_size Font size of the row names (default: 7)
#' @param column_names_size Font size of the column names (default: 9)
#' @param column_names_rotation Rotation of the column names (default: 45)
#' @param column_dend_height Height of the column dendrogram (a number in mm, or a unit)
#' @param row_dend_width Width of the row dendrogram (a number in mm, or a unit)
#' @param show_heatmap_legend Show the heatmap legend (default: TRUE)
#' @param show_annotation_legend Show the annotation legends (default: TRUE)
#' @param border_color Border colour of the cells (NULL, TRUE, or a colour)
#' @param heatmap_title Main heatmap title (default: NULL, no title).
#'   "\{mode\}" can be used as a placeholder, and is replaced by the mode name
#' @param heatmap_title_size Font size of the main title (default: 14)
#' @param heatmap_title_face Font face of the title (default: "bold")
#' @param export_path Base path used to export the data to TSV (default: NULL).
#'   The mode name is appended to the file name (e.g. "export_all.tsv", "export_B-A.tsv")
#' @param export_modes Vector of modes to export (default: NULL, exports all of them).
#'   Only relevant when export_path is set. Example: c("all", "B-A")
#' @param export_file Base path used to export the plots. The mode is appended to the name.
#'   The format is taken from the extension (.png, .svg, .pdf)
#' @param plot_width Plot width in inches (default: 10)
#' @param plot_height Plot height in inches (default: 8)
#' @param export_dpi Resolution for PNG output in dots per inch (default: 300)
#'
#' @return Named list of tidyHeatmap objects
#'
#' @examples
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("tidyHeatmap", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#'   # "any" plus one heatmap per named comparison
#'   hm_list <- proteomics_heatmap_list(
#'     res$PCA_Input,
#'     modes           = c("any", "B-A"),
#'     sample_order    = "condition",
#'     condition_order = c("A", "B", "D"),
#'     palette_value   = "brewer:RdYlBu"
#'   )
#'   print(names(hm_list))   # print(hm_list[["B-A"]]) draws one of them
#' }
#'
#' @export
proteomics_heatmap_list <- function(data,
                                    modes = c("all", "any"),
                                    alpha = 0.05,
                                    feature_ids = NULL,
                                    row_annotation = NULL,
                                    row_annotation_cols = NULL,
                                    row_annotation_palette = NULL,
                                    row_annotation_size = NULL,
                                    row_annotation_name_size = 8,
                                    row_order_by = NULL,
                                    split_rows_by = NULL,
                                    scale_data = c("row", "none", "column"),
                                    sample_order = "clustering",
                                    condition_order = NULL,
                                    cluster_rows = FALSE,
                                    cluster_columns = TRUE,
                                    show_row_names = NULL,
                                    show_column_names = TRUE,
                                    palette_value = NULL,
                                    palette_annotation = NULL,
                                    reverse_palette = FALSE,
                                    row_title_size = 10,
                                    column_title_size = 10,
                                    show_row_title = TRUE,
                                    show_column_title = TRUE,
                                    show_annotation = TRUE,
                                    split_by_condition = FALSE,
                                    show_adjp_annotation = TRUE,
                                    palette_adjp = NULL,
                                    row_names_size = 7,
                                    column_names_size = 9,
                                    column_names_rotation = 45,
                                    column_dend_height = NULL,
                                    row_dend_width = NULL,
                                    show_heatmap_legend = TRUE,
                                    show_annotation_legend = TRUE,
                                    border_color = NULL,
                                    heatmap_title = NULL,
                                    heatmap_title_size = 14,
                                    heatmap_title_face = "bold",
                                    export_path = NULL,
                                    export_modes = NULL,
                                    export_file = NULL,
                                    plot_width = 10,
                                    plot_height = 8,
                                    export_dpi = 300) {

  scale_data <- match.arg(scale_data)

  # ---------------------------------------------------------------------------
  # 1) Input validation
  # ---------------------------------------------------------------------------

  required_cols <- c("SampleID", "FeatureID", "Intensity", "Condition", "Replicate")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # 2) Detect the available comparisons (adjP_* columns)
  # ---------------------------------------------------------------------------

  adjp_cols <- grep("^adjP_", names(data), value = TRUE)
  available_comparisons <- sub("^adjP_", "", adjp_cols)

  # ---------------------------------------------------------------------------
  # 3) Build a heatmap for each mode
  # ---------------------------------------------------------------------------

  hm_list <- list()

  for (m in modes) {

    # Resolve the internal mode and the comparison (where applicable)
    if (m == "all") {
      internal_mode <- "all"
      comparison <- NULL
      row_title <- "All Proteins"
    } else if (m == "any") {
      internal_mode <- "any"
      comparison <- NULL
      row_title <- "DEPs (any comparison)"
    } else {
      # It is a specific comparison (mode "target")
      if (!(m %in% available_comparisons)) {
        warning("Comparison '", m, "' not found. It will be skipped.")
        next
      }
      internal_mode <- "target"
      comparison <- m
      row_title <- paste0("DEPs (", m, ")")
    }

    # Resolve the title (substitute {mode} when present)
    current_title <- NULL
    if (!is.null(heatmap_title)) {
      current_title <- gsub("\\{mode\\}", m, heatmap_title)
    }

    # Resolve export_path (append the mode to the file name)
    current_export_path <- NULL
    if (!is.null(export_path) && nzchar(export_path)) {
      # Check whether this mode should be exported
      should_export <- is.null(export_modes) || m %in% export_modes
      if (should_export) {
        # Split off the directory, the base name and the extension
        dir_path <- dirname(export_path)
        base_name <- tools::file_path_sans_ext(basename(export_path))
        ext <- tools::file_ext(export_path)
        if (nzchar(ext)) ext <- paste0(".", ext) else ext <- ".tsv"
        current_export_path <- file.path(dir_path, paste0(base_name, "_", m, ext))
      }
    }

    # Resolve export_file (append the mode to the plot file name)
    current_export_file <- NULL
    if (!is.null(export_file) && nzchar(export_file)) {
      # Check whether this mode should be exported
      should_export <- is.null(export_modes) || m %in% export_modes
      if (should_export) {
        # Split off the directory, the base name and the extension
        dir_path <- dirname(export_file)
        base_name <- tools::file_path_sans_ext(basename(export_file))
        ext <- tools::file_ext(export_file)
        if (nzchar(ext)) ext <- paste0(".", ext) else ext <- ".png"
        current_export_file <- file.path(dir_path, paste0(base_name, "_", m, ext))
      }
    }

    # Try to build the heatmap (it can fail when there are not enough proteins)
    hm <- tryCatch({
      proteomics_heatmap(
        data = data,
        mode = internal_mode,
        alpha = alpha,
        comparison = comparison,
        feature_ids = feature_ids,
        row_annotation = row_annotation,
        row_annotation_cols = row_annotation_cols,
        row_annotation_palette = row_annotation_palette,
        row_annotation_size = row_annotation_size,
        row_annotation_name_size = row_annotation_name_size,
        row_order_by = row_order_by,
        split_rows_by = split_rows_by,
        scale_data = scale_data,
        sample_order = sample_order,
        condition_order = condition_order,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        palette_value = palette_value,
        palette_annotation = palette_annotation,
        reverse_palette = reverse_palette,
        row_title = row_title,
        column_title = "Samples",
        row_title_size = row_title_size,
        column_title_size = column_title_size,
        show_row_title = show_row_title,
        show_column_title = show_column_title,
        show_annotation = show_annotation,
        split_by_condition = split_by_condition,
        show_adjp_annotation = show_adjp_annotation,
        palette_adjp = palette_adjp,
        row_names_size = row_names_size,
        column_names_size = column_names_size,
        column_names_rotation = column_names_rotation,
        column_dend_height = column_dend_height,
        row_dend_width = row_dend_width,
        show_heatmap_legend = show_heatmap_legend,
        show_annotation_legend = show_annotation_legend,
        border_color = border_color,
        heatmap_title = current_title,
        heatmap_title_size = heatmap_title_size,
        heatmap_title_face = heatmap_title_face,
        export_path = current_export_path,
        export_file = current_export_file,
        plot_width = plot_width,
        plot_height = plot_height,
        export_dpi = export_dpi
      )
    }, error = function(e) {
      warning("Could not build the heatmap for '", m, "': ", e$message)
      return(NULL)
    })

    if (!is.null(hm)) {
      hm_list[[m]] <- hm
    }
  }

  hm_list
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# --- Load the data ---
# hm_input <- arrow::read_parquet("PCA_Input.parquet")
# hm_input <- readr::read_tsv("PCA_Input.tsv")

# --- Basic example: heatmap of every protein ---
# hm_all <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_all

# --- Heatmap with column clustering ---
# hm_cluster <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = "clustering",
#   cluster_columns = TRUE,
#   cluster_rows = FALSE
# )
# hm_cluster

# --- Heatmap of the proteins significant in any comparison ---
# hm_any <- proteomics_heatmap(
#   data = hm_input,
#   mode = "any",
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_any

# --- Heatmap of the proteins significant in one specific comparison ---
# hm_target <- proteomics_heatmap(
#   data = hm_input,
#   mode = "target",
#   comparison = "B-A",
#   alpha = 0.05,
#   scale_data = "row"
# )
# hm_target

# --- With a custom paletteer palette ---
# hm_viridis <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   palette_value = "viridis::viridis",
#   reverse_palette = TRUE
# )
# hm_viridis

# --- With an RColorBrewer palette ---
# hm_brewer <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   palette_value = "brewer:RdYlBu",
#   palette_annotation = "brewer:Set1"
# )
# hm_brewer

# --- Build several heatmaps at once ---
# hm_list <- proteomics_heatmap_list(
#   data = hm_input,
#   modes = c("all", "any", "B-A", "C-A", "D-A"),
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_list[["all"]]
# hm_list[["any"]]
# hm_list[["B-A"]]

# --- Custom sample order ---
# custom_order <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2")
# hm_custom <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = custom_order,
#   cluster_columns = FALSE
# )
# hm_custom

# --- Heatmap with a custom title ---
# hm_title <- proteomics_heatmap(
#   data = hm_input,
#   mode = "any",
#   scale_data = "row",
#   heatmap_title = "Differential Expression Analysis",
#   heatmap_title_size = 16,
#   heatmap_title_face = "bold"
# )
# hm_title

# --- List of heatmaps with dynamic titles ---
# hm_list <- proteomics_heatmap_list(
#   data = hm_input,
#   modes = c("all", "any", "B-A"),
#   heatmap_title = "Proteomics Heatmap: {mode}",
#   heatmap_title_size = 14
# )
# # The titles will be: "Proteomics Heatmap: all", "Proteomics Heatmap: any", "Proteomics Heatmap: B-A"

# --- Available palettes ---
# Diverging (good for scaled data):
#   - NULL (default RdBu-like)
#   - "brewer:RdBu", "brewer:RdYlBu", "brewer:PiYG", "brewer:BrBG"
#   - "viridis::plasma", "viridis::inferno"
#
# Sequential (good for unscaled data):
#   - "viridis::viridis", "viridis::magma"
#   - "brewer:Blues", "brewer:Reds", "brewer:YlOrRd"
#
# For annotations:
#   - "brewer:Set1", "brewer:Set2", "brewer:Dark2"
#   - "ggsci::category10_d3", "ggsci::nrc_npg"
