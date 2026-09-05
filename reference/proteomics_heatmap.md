# Create a tidyHeatmap Heatmap for Proteomics

Create a tidyHeatmap Heatmap for Proteomics

## Usage

``` r
proteomics_heatmap(
  data,
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
  export_dpi = 300
)
```

## Arguments

- data:

  Data frame in long format (see prepare_heatmap_data for the structure)

- mode:

  Protein filtering mode: "all", "any", or "target"

- alpha:

  Significance threshold for DEP proteins (default: 0.05)

- comparison:

  Name of the comparison for mode "target" (e.g. "B-A")

- feature_ids:

  Vector of specific FeatureIDs to show (default: NULL, uses the
  mode-based filtering). When supplied, only these proteins are shown
  and the mode/alpha filtering is ignored.

- row_annotation:

  Row annotations (FeatureIDs). It can be:

  - Path to a TSV file with a "FeatureID" column plus additional
    categorical columns

  - Data frame with the same structure

  - NULL: no row annotations (default)

- row_annotation_cols:

  Vector of column names to display as annotations. Default: NULL (uses
  every column except FeatureID)

- row_annotation_palette:

  Named list of palettes, one per annotation. Example: list(Pathway =
  "brewer:Set1", Function = c("red", "blue", "green"))

- row_annotation_size:

  Width of the row annotation bars. It can be:

  - Number: interpreted as centimetres (e.g. 0.3 = 0.3cm)

  - A unit object: grid::unit(0.3, "cm")

  - NULL: use tidyHeatmap's default value

- row_annotation_name_size:

  Font size of the row annotation names (default: 8)

- row_order_by:

  Row order. It can be:

  - NULL: no particular ordering (default)

  - "clustering": order by hierarchical clustering

  - A column name: order by that annotation (e.g. "Specie")

  - A vector of columns: order sequentially (e.g. c("Specie", "Process",
    "Function"))

  - Include "adjP" when mode = "target" to order by adjusted p-value

- split_rows_by:

  Name of the annotation column used to split the rows into groups
  (default: NULL)

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

  Show the row names (default: TRUE when \<= 50 proteins)

- show_column_names:

  Show the column names (default: TRUE)

- palette_value:

  Palette for the heatmap values (see get_heatmap_palette)

- palette_annotation:

  Palette for the Condition annotation

- reverse_palette:

  Reverse the value palette (default: FALSE)

- row_title:

  Title for the rows

- column_title:

  Title for the columns

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

  Show the adjP row annotation when mode="target" (default: TRUE)

- palette_adjp:

  Palette for the adjP annotation. It can be:

  - NULL: use the default palette (red-orange-white)

  - A vector of 3 colours: c(color_0, color_middle, color_0.05)

  - String "brewer:name": use an RColorBrewer palette

  - String "package::palette": use a paletteer palette

- row_names_size:

  Font size of the row names (default: 7)

- column_names_size:

  Font size of the column names (default: 9)

- column_names_rotation:

  Rotation of the column names in degrees (default: 45)

- column_dend_height:

  Height of the column dendrogram. It can be:

  - Number: interpreted as millimetres (e.g. 30 = 30mm)

  - A unit object: grid::unit(2, "cm")

  - NULL: use ComplexHeatmap's default value

- row_dend_width:

  Width of the row dendrogram. Same format as column_dend_height

- show_heatmap_legend:

  Show the heatmap legend (default: TRUE)

- show_annotation_legend:

  Show the annotation legends (default: TRUE)

- border_color:

  Border colour of the heatmap cells. It can be:

  - NULL or FALSE: no border (default)

  - TRUE: black border

  - A colour string: that specific colour (e.g. "black", "grey",
    "#CCCCCC")

- heatmap_title:

  Main heatmap title (default: NULL, no title)

- heatmap_title_size:

  Font size of the main title (default: 14)

- heatmap_title_face:

  Font face of the title: "plain", "bold", "italic", "bold.italic"
  (default: "bold")

- export_path:

  Path used to export the heatmap data to TSV (default: NULL, nothing is
  exported). The file contains FeatureID, the per-sample intensity
  values, and the metadata (adjP where applicable).

- export_file:

  Path used to export the plot. The format is taken from the extension:

  - .png: PNG image (raster)

  - .svg: SVG image (vector)

  - .pdf: PDF document (vector)

- plot_width:

  Plot width in inches (default: 10). To get pixels: pixels = inches x
  dpi (e.g. 10" x 300dpi = 3000px)

- plot_height:

  Plot height in inches (default: 8). To get pixels: pixels = inches x
  dpi (e.g. 8" x 300dpi = 2400px)

- export_dpi:

  Resolution for PNG output in dots per inch (default: 300). Higher dpi
  = more detail. Common values: 72 (web), 150 (draft), 300 (publication)

## Value

tidyHeatmap/ComplexHeatmap object

## Examples

``` r
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("tidyHeatmap", quietly = TRUE)) {
  data(nadia_dia)
  res <- process_proteomics(nadia_dia, verbose = FALSE)

  # Proteins significant in the B-A comparison, samples grouped by condition
  hm <- proteomics_heatmap(
    res$PCA_Input,
    mode            = "target",
    comparison      = "B-A",
    scale_data      = "row",
    sample_order    = "condition",
    condition_order = c("A", "B", "D")
  )
  print(class(hm))   # print(hm) itself draws it on the current device
}
#> [1] "InputHeatmap"
#> attr(,"package")
#> [1] "tidyHeatmap"
```
