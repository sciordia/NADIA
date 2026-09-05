# Interactive Highcharts Boxplot for Proteomics

Interactive Highcharts Boxplot for Proteomics

## Usage

``` r
boxplot_highchart_list(
  data,
  assays = NULL,
  color_by = "Condition",
  group_order = NULL,
  palette = NULL,
  title = NULL,
  subtitle = NULL,
  show_outliers = TRUE,
  outlier_jitter = 0.15,
  outlier_size = 3,
  box_width = 20,
  horizontal = TRUE,
  height = NULL
)
```

## Arguments

- data:

  Data frame with columns: Column, Assay, Intensity, Condition

- assays:

  Vector of assays to include (NULL = all of them)

- color_by:

  Column used for colouring (default: "Condition")

- group_order:

  Order of the groups/conditions

- palette:

  Colour palette: "ggsci::palette", "brewer:Name", or a vector

- title:

  Chart title (optional). Use `{assay}` as a placeholder (e.g. "Boxplot:
  {assay}" -\> "Boxplot: ImpSeqRob_Min")

- subtitle:

  Chart subtitle (optional). Use `{assay}` as a placeholder

- show_outliers:

  Show outliers lying outside the whiskers (default: TRUE)

- outlier_jitter:

  Amount of horizontal jitter applied to the outliers (default: 0.15)

- outlier_size:

  Radius of the outlier points (default: 3)

- box_width:

  Width of the boxplot boxes in pixels (default: 20)

- horizontal:

  Horizontal orientation (default: TRUE)

- height:

  Chart height in pixels

## Value

List of highchart objects (one per assay)

## Examples

``` r
data(nadia_dia)
res <- process_proteomics(nadia_dia, verbose = FALSE)

# One boxplot per assay; keep only the imputed one
hc_list <- boxplot_highchart_list(
  res$BoxPlot_Input,
  assays      = "Impseqrob_min",
  group_order = c("A", "B", "D"),
  title       = "Boxplot: {assay}"
)
names(hc_list)
#> [1] "Impseqrob_min"
```
