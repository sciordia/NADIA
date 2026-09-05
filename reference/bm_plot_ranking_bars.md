# Bar chart of final ranks

Horizontal bars of rank_final by Assay, ordered best to worst.

## Usage

``` r
bm_plot_ranking_bars(ranking_result, type = "mean", title = NULL)
```

## Arguments

- ranking_result:

  List returned by
  [`bm_compute_ranking()`](https://sciordia.github.io/NADIA/reference/bm_compute_ranking.md)

- type:

  "mean" or "median" aggregation (default: "mean")

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
ranking <- bm_compute_ranking(opdea_all)
bm_plot_ranking_bars(ranking)
```
