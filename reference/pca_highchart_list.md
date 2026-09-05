# Build a List of PCA Plots for Multiple Subsets

Automatically builds PCA plots for "all", "any" and/or specific
comparisons. Similar to factoextra::fviz_pca_ind with ellipse options.

## Usage

``` r
pca_highchart_list(
  pca_input,
  modes = c("all", "any"),
  alpha = 0.05,
  color_by = "Condition",
  group_order = NULL,
  palette = NULL,
  addEllipses = TRUE,
  ellipse_type = c("convex", "confidence"),
  ellipse_level = 0.95,
  ellipse_fill_opacity = 0.12,
  ellipse_line_width = 1,
  ellipse_npoints = 100,
  point_size = 5,
  show_labels = FALSE,
  label_size = 10,
  center = TRUE,
  scale. = TRUE,
  filter_samples_to_comparison = FALSE
)
```

## Arguments

- pca_input:

  Data frame in long format (see build_pca_scores for the structure)

- modes:

  Vector of modes to generate: "all", "any", and/or comparison names
  (default: c("all", "any"))

- alpha:

  Significance threshold for DEPs (default: 0.05)

- color_by:

  Column used to colour the points (default: "Condition")

- group_order:

  Order of the groups/conditions (optional)

- palette:

  Named vector of colours, or NULL for an automatic palette

- addEllipses:

  Show ellipses/hulls around the groups (default: TRUE)

- ellipse_type:

  Ellipse type: "convex" or "confidence" (default: "convex")

- ellipse_level:

  Confidence level for ellipse_type = "confidence" (default: 0.95)

- ellipse_fill_opacity:

  Ellipse fill opacity (0-1, default: 0.12)

- ellipse_line_width:

  Outline line width (default: 1)

- ellipse_npoints:

  Number of points for the confidence ellipse (default: 100)

- point_size:

  Point radius (default: 5)

- show_labels:

  Show the point labels (SampleID) (default: FALSE)

- label_size:

  Label font size in px (default: 10)

- center:

  Center the data before PCA (default: TRUE)

- scale.:

  Scale the data before PCA (default: TRUE)

- filter_samples_to_comparison:

  For specific comparisons, restrict the samples to the conditions
  involved (default: FALSE)

## Value

Named list of highchart objects

## Examples

``` r
data(nadia_dia)
res <- process_proteomics(nadia_dia, verbose = FALSE)

# One PCA over all proteins and one over the significant ones
hc_pcas <- pca_highchart_list(
  res$PCA_Input,
  modes         = c("all", "any"),
  group_order   = c("A", "B", "D"),
  ellipse_type  = "confidence",
  ellipse_level = 0.95
)
names(hc_pcas)
#> [1] "all" "any"
```
