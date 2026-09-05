# Build a List of Heatmaps for Several Subsets

Automatically builds heatmaps for "all", "any" and/or specific
comparisons.

## Usage

``` r
proteomics_heatmap_list(
  data,
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
  export_dpi = 300
)
```

## Arguments

- data:

  Data frame in long format (see prepare_heatmap_data for the structure)

- modes:

  Vector of modes to generate: "all", "any", and/or comparison names
  (default: c("all", "any"))

- alpha:

  Significance threshold for DEP proteins (default: 0.05)

- feature_ids:

  Vector of specific FeatureIDs to show (default: NULL)

- row_annotation:

  Row annotations (TSV path or data.frame)

- row_annotation_cols:

  Columns to use as row annotations

- row_annotation_palette:

  List of palettes for the row annotations

- row_annotation_size:

  Width of the row annotation bars (a number in cm, or a unit)

- row_annotation_name_size:

  Font size of the row annotation names (default: 8)

- row_order_by:

  Row order: "clustering", annotation column(s), or include "adjP" for
  mode="target"

- split_rows_by:

  Annotation column used to split the rows into groups

- scale_data:

  Scaling type: "none", "row", "column" (default: "row")

- sample_order:

  Sample order: "clustering", "condition", or a custom vector

- condition_order:

  Condition order when sample_order = "condition"

- cluster_rows:

  Cluster the rows (default: FALSE)

- cluster_columns:

  Cluster the columns (default: TRUE)

- show_row_names:

  Show the row names (default: auto)

- show_column_names:

  Show the column names (default: TRUE)

- palette_value:

  Palette for the heatmap values

- palette_annotation:

  Palette for the Condition annotation

- reverse_palette:

  Reverse the value palette (default: FALSE)

- row_title_size:

  Font size of the row title (default: 10)

- column_title_size:

  Font size of the column title (default: 10)

- show_row_title:

  Show the row title (default: TRUE)

- show_column_title:

  Show the column title (default: TRUE)

- show_annotation:

  Show the Condition annotation (default: TRUE)

- split_by_condition:

  Split the heatmap by condition (default: FALSE)

- show_adjp_annotation:

  Show the adjP row annotation for comparisons (default: TRUE)

- palette_adjp:

  Palette for the adjP annotation (see proteomics_heatmap)

- row_names_size:

  Font size of the row names (default: 7)

- column_names_size:

  Font size of the column names (default: 9)

- column_names_rotation:

  Rotation of the column names (default: 45)

- column_dend_height:

  Height of the column dendrogram (a number in mm, or a unit)

- row_dend_width:

  Width of the row dendrogram (a number in mm, or a unit)

- show_heatmap_legend:

  Show the heatmap legend (default: TRUE)

- show_annotation_legend:

  Show the annotation legends (default: TRUE)

- border_color:

  Border colour of the cells (NULL, TRUE, or a colour)

- heatmap_title:

  Main heatmap title (default: NULL, no title). "{mode}" can be used as
  a placeholder, and is replaced by the mode name

- heatmap_title_size:

  Font size of the main title (default: 14)

- heatmap_title_face:

  Font face of the title (default: "bold")

- export_path:

  Base path used to export the data to TSV (default: NULL). The mode
  name is appended to the file name (e.g. "export_all.tsv",
  "export_B-A.tsv")

- export_modes:

  Vector of modes to export (default: NULL, exports all of them). Only
  relevant when export_path is set. Example: c("all", "B-A")

- export_file:

  Base path used to export the plots. The mode is appended to the name.
  The format is taken from the extension (.png, .svg, .pdf)

- plot_width:

  Plot width in inches (default: 10)

- plot_height:

  Plot height in inches (default: 8)

- export_dpi:

  Resolution for PNG output in dots per inch (default: 300)

## Value

Named list of tidyHeatmap objects

## Examples

``` r
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("tidyHeatmap", quietly = TRUE)) {
  data(nadia_dia)
  res <- process_proteomics(nadia_dia, verbose = FALSE)

  # "any" plus one heatmap per named comparison
  hm_list <- proteomics_heatmap_list(
    res$PCA_Input,
    modes           = c("any", "B-A"),
    sample_order    = "condition",
    condition_order = c("A", "B", "D"),
    palette_value   = "brewer:RdYlBu"
  )
  print(names(hm_list))   # print(hm_list[["B-A"]]) draws one of them
}
#> [1] "any" "B-A"
```
