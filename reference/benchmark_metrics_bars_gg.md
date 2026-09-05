# Grouped Metrics Bar Chart (ggplot2)

Grouped bar chart: Comparisons x multiple metrics side by side.

## Usage

``` r
benchmark_metrics_bars_gg(
  metrics_table,
  metrics_to_show = c("Sensitivity", "Specificity", "Precision", "F1", "Accuracy"),
  title = "Benchmark Metrics by Comparison",
  bar_width = 0.7
)
```

## Arguments

- metrics_table:

  Data frame from compute_benchmark_metrics()

- metrics_to_show:

  Character vector of metrics to display

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
metrics <- compute_benchmark_metrics(de, ev)
benchmark_metrics_bars_gg(metrics, metrics_to_show = c("Sensitivity", "F1"))
```
