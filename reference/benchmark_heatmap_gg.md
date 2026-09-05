# Performance Heatmap (ggplot2)

Heatmap of Comparisons x Metrics (AUC, Sensitivity, Specificity,
Precision, F1, Accuracy). Gradient red -\> yellow -\> green.

## Usage

``` r
benchmark_heatmap_gg(
  metrics_table,
  metrics_to_show = c("AUC", "Sensitivity", "Specificity", "Precision", "F1", "Accuracy"),
  title = "Benchmark Performance Heatmap",
  text_size = 4,
  axis_text_size = 11
)
```

## Arguments

- metrics_table:

  Data frame from compute_benchmark_metrics()

- metrics_to_show:

  Character vector of metrics to display

- title:

  Plot title

- text_size:

  Size of cell text labels

- axis_text_size:

  Size of axis text

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
benchmark_heatmap_gg(metrics, metrics_to_show = c("AUC", "F1", "MCC"))
```
