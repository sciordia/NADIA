# AUC Bar Chart (ggplot2)

Horizontal bar chart of AUC per comparison, colored by value, with a
reference line at 0.5.

## Usage

``` r
benchmark_auc_bars_gg(
  metrics_table,
  title = "AUC by Comparison",
  bar_width = 0.7
)
```

## Arguments

- metrics_table:

  Data frame from compute_benchmark_metrics()

- title:

  Plot title

- bar_width:

  Bar width (default: 0.7)

## Value

ggplot2 object

## Examples

``` r
data(nadia_dia)
de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results
sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                    package = "NADIA"))
de$Species <- sp$PG.OrganismId[match(de$Protein.IDs, sp$PG.ProteinGroups)]
ev <- data.frame(Comparison = rep(c("B-A", "D-A"), each = 2),
                 Species = c("ECOLI", "YEAST"),
                 expected_logFC = c(1, -0.58, 2, -3.3))
if (requireNamespace("pROC", quietly = TRUE)) {
  metrics <- compute_benchmark_metrics(de, ev)
  benchmark_auc_bars_gg(metrics)
}
```
