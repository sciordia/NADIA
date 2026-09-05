# Stacked bar chart of TP/FP/FN/TN (OpDEA Figure 5 style)

Horizontal stacked bars showing confusion matrix counts per method. If
`comparison` is NULL, generates a faceted plot with all comparisons.

## Usage

``` r
bm_plot_confusion_stacked(confusion_combined, comparison = NULL, title = NULL)
```

## Arguments

- confusion_combined:

  data.frame with columns: Assay, Comparison, TP, FP, TN, FN

- comparison:

  Optional: single comparison to plot. If NULL, faceted plot with all
  comparisons.

- title:

  Optional plot title

## Value

ggplot2 object

## Examples

``` r
data(nadia_dia)
sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                    package = "NADIA"))
ev <- data.frame(Comparison = rep(c("B-A", "D-A"), each = 2),
                 Species = c("ECOLI", "YEAST"),
                 expected_logFC = c(1, -0.58, 2, -3.3))
bench_all <- do.call(rbind, lapply(c("cycloess", "quantile"), function(nm) {
  de <- process_proteomics(nadia_dia, norm_method = nm,
                           verbose = FALSE)$DEPs_results
  de$Species <- sp$PG.OrganismId[match(de$Protein.IDs, sp$PG.ProteinGroups)]
  cbind(compute_benchmark_metrics(de, ev), Assay = nm)
}))
bm_plot_confusion_stacked(bench_all)

bm_plot_confusion_stacked(bench_all, comparison = "B-A")
```
