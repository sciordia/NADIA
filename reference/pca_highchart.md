# Interactive PCA Plot with Highcharts

Builds a PCA scatter plot with the option of showing per-group ellipses
or convex hulls, similar to factoextra::fviz_pca_ind.

## Usage

``` r
pca_highchart(
  scores_df,
  color_by = "Condition",
  group_order = NULL,
  palette = NULL,
  title = NULL,
  addEllipses = TRUE,
  ellipse_type = c("convex", "confidence"),
  ellipse_level = 0.95,
  ellipse_fill_opacity = 0.12,
  ellipse_line_width = 1,
  ellipse_npoints = 100,
  point_size = 5,
  show_labels = FALSE,
  label_size = 10
)
```

## Arguments

- scores_df:

  Data frame produced by build_pca_scores() with columns: PC1, PC2,
  PC1_Perc, PC2_Perc, Subset, Condition, Replicate, SampleID

- color_by:

  Column used to colour the points (default: "Condition")

- group_order:

  Vector with the order of the groups/conditions (optional)

- palette:

  Named vector of colours, or NULL for an automatic palette

- title:

  Chart title (optional, defaults to Subset)

- addEllipses:

  Show ellipses/hulls around the groups (default: TRUE)

- ellipse_type:

  Ellipse type: "convex" for a convex hull or "confidence" for a
  confidence ellipse based on the normal distribution (default:
  "convex")

- ellipse_level:

  Confidence level for ellipse_type = "confidence" (default: 0.95).
  Typical values: 0.95, 0.90, 0.68

- ellipse_fill_opacity:

  Fill opacity (0-1, default: 0.12)

- ellipse_line_width:

  Outline line width (default: 1)

- ellipse_npoints:

  Number of points used to draw the confidence ellipse (default: 100).
  Only applies to ellipse_type = "confidence"

- point_size:

  Point radius (default: 5)

- show_labels:

  Show the point labels (SampleID) (default: FALSE)

- label_size:

  Label font size in px (default: 10)

## Value

A highchart object

## Examples

``` r
# The scores are normally built internally by pca_highchart_list();
# this is the layout the function expects (one row per sample).
scores <- data.frame(
  SampleID  = paste(rep(c("A", "B", "D"), each = 4), 1:4, sep = "_"),
  Condition = rep(c("A", "B", "D"), each = 4),
  Replicate = rep(1:4, times = 3),
  PC1 = c(-5.1, -4.6, -5.4, -4.9, 0.4, 0.9, 0.2, 0.7, 4.8, 5.2, 4.5, 5.0),
  PC2 = c(1.2, 0.7, -0.9, -1.1, 2.4, 1.8, -1.6, -2.2, 0.9, 0.4, -0.8, -1.3),
  PC1_Perc = 62.4,
  PC2_Perc = 18.7,
  Subset = "All proteins",
  stringsAsFactors = FALSE
)

hc <- pca_highchart(scores, group_order = c("A", "B", "D"),
                    ellipse_type = "confidence", show_labels = TRUE)
class(hc)
#> [1] "highchart"  "htmlwidget"
```
