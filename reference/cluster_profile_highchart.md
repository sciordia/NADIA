# Cluster profile plot with Highcharts

Builds an interactive plot showing the expression profiles of the
proteins in a given cluster.

## Usage

``` r
cluster_profile_highchart(
  data,
  cluster,
  conditions = NULL,
  min_membership = NULL,
  show_centroid = TRUE,
  centroid_summary = c("mean", "median"),
  cluster_color = NULL,
  line_width = 1,
  line_opacity = 0.4,
  centroid_width = 3,
  title = NULL,
  height = NULL
)
```

## Arguments

- data:

  Pattern Profiler DataFrame (LONG format)

- cluster:

  Number of the cluster to plot

- conditions:

  Vector of condition names (order for the X axis)

- min_membership:

  Additional membership filter (NULL = no filter)

- show_centroid:

  Show the centroid line (default: TRUE)

- centroid_summary:

  Method for the centroid: "mean" or "median"

- cluster_color:

  Cluster colour (NULL = automatic)

- line_width:

  Width of the profile lines (default: 1)

- line_opacity:

  Opacity of the lines (default: 0.4)

- centroid_width:

  Width of the centroid line (default: 3)

- title:

  Custom title (optional)

- height:

  Plot height in pixels

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

  hc <- cluster_profile_highchart(pp$long_output, cluster = 1,
                                  conditions = c("A", "B", "D"))
  print(class(hc))
}
#> [1] "highchart"  "htmlwidget"
```
