# Boxplot comparison of metric distributions across methods

Faceted boxplot showing the distribution of each metric across
comparisons for every method.

## Usage

``` r
bm_plot_metrics_comparison(opdea_combined, metrics = .BM_METRICS, title = NULL)
```

## Arguments

- opdea_combined:

  data.frame with Assay, Comparison, and metric columns

- metrics:

  Character vector of metrics to plot (default: all 5)

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
opdea_all <- do.call(rbind, lapply(c("cycloess", "quantile"), function(nm) {
  de <- process_proteomics(nadia_dia, norm_method = nm,
                           verbose = FALSE)$DEPs_results
  de$Species <- sp$PG.OrganismId[match(de$Protein.IDs, sp$PG.ProteinGroups)]
  cbind(compute_opdea_metrics(de, ev), Assay = nm)
}))
bm_plot_metrics_comparison(opdea_all)
```
