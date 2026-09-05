# Bar chart of final ranks for a single comparison

Bar chart of final ranks for a single comparison

## Usage

``` r
bm_plot_ranking_bars_by_comp(ranking_by_comp, comparison, title = NULL)
```

## Arguments

- ranking_by_comp:

  List returned by
  [`bm_compute_ranking_by_comparison()`](https://sciordia.github.io/NADIA/reference/bm_compute_ranking_by_comparison.md)

- comparison:

  Character string: which comparison to plot

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
by_comp <- bm_compute_ranking_by_comparison(opdea_all)
bm_plot_ranking_bars_by_comp(by_comp, comparison = "D-A")
```
