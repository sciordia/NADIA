# Interactive Volcano Plot with Highcharter

Interactive Volcano Plot with Highcharter

## Usage

``` r
volcano_highchart_list(
  de_res,
  ain = NULL,
  comparisons = NULL,
  lfc_thr = 0,
  alpha = 0.05,
  p_col = "adj.P.Val",
  point_size = 4,
  colors = NULL,
  palette = NULL,
  show_top_genes = 0,
  highlight_genes = NULL,
  title = NULL
)
```

## Arguments

- de_res:

  Data frame with differential expression results

- ain:

  Vector of assays to filter on (optional)

- comparisons:

  Vector of comparisons to include (optional)

- lfc_thr:

  log2 fold-change threshold (default: 0)

- alpha:

  Significance threshold (default: 0.05)

- p_col:

  Column of p-values to use

- point_size:

  Point size (default: 4)

- colors:

  List of colours for "up", "down", "ns" (ignored when palette is used)

- palette:

  Name of a paletteer palette (e.g. "ggsci::default_jco"). Uses 3
  colours: up, down, ns

- show_top_genes:

  Number of top genes to label by significance (default: 0)

- highlight_genes:

  Vector of gene names to highlight manually (default: NULL)

- title:

  Custom chart title (default: NULL, uses "Comparison (Assay)"). Use
  `{comparison}` as a placeholder (e.g. "Volcano Plot: {comparison}" -\>
  "Volcano Plot: B-A")

## Value

List of highchart objects

## Examples

``` r
data(nadia_dia)
res <- process_proteomics(nadia_dia, verbose = FALSE)

# One volcano per comparison, labelling the 10 most significant proteins
hc_list <- volcano_highchart_list(
  res$DEPs_results,
  comparisons     = c("B-A", "D-A"),
  lfc_thr         = 1,
  alpha           = 0.05,
  show_top_genes  = 10
)
names(hc_list)
#> [1] "B-A" "D-A"
```
