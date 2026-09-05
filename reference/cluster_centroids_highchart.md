# Centroid plot for all the clusters

Builds a comparative plot with the centroids of all the clusters.

## Usage

``` r
cluster_centroids_highchart(
  data,
  conditions = NULL,
  clusters = NULL,
  min_membership = NULL,
  centroid_summary = c("mean", "median"),
  palette = NULL,
  line_width = 2.5,
  show_markers = TRUE,
  title = NULL,
  height = NULL
)
```

## Arguments

- data:

  Pattern Profiler DataFrame (LONG format)

- conditions:

  Vector of condition names

- clusters:

  Clusters to include (NULL = all)

- min_membership:

  Membership filter applied before computing the centroids

- centroid_summary:

  Method: "mean" or "median"

- palette:

  Colour palette

- line_width:

  Line width (default: 2.5)

- show_markers:

  Show markers on the points (default: TRUE)

- title:

  Custom title

- height:

  Plot height

## Value

highchart object

## Examples

``` r
if (requireNamespace("Mfuzz", quietly = TRUE) &&
    requireNamespace("e1071", quietly = TRUE)) {
  data(nadia_dia)
  res <- process_proteomics(nadia_dia, verbose = FALSE)
  pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
                                  assay_name = "Impseqrob_min",
                                  auto_select_c = FALSE, c = 3,
                                  verbose = FALSE)

  # All the cluster centroids on one chart
  hc <- cluster_centroids_highchart(pp$long_output,
                                    conditions = c("A", "B", "D"),
                                    centroid_summary = "median")
  print(class(hc))
}
#> Warning: no DISPLAY variable so Tk is not available
#> [1] "highchart"  "htmlwidget"
```
