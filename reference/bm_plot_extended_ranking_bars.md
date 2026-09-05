# Horizontal bar chart of extended combined ranking

Horizontal bar chart of extended combined ranking

## Usage

``` r
bm_plot_extended_ranking_bars(extended_ranking, title = NULL)
```

## Arguments

- extended_ranking:

  List from bm_compute_extended_ranking()

- title:

  Optional plot title

## Value

ggplot object

## Examples

``` r
data(nadia_dia)
sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                    package = "NADIA"))
ev <- data.frame(Comparison = rep(c("B-A", "D-A"), each = 2),
                 Species = c("ECOLI", "YEAST"),
                 expected_logFC = c(1, -0.58, 2, -3.3))
opdea_all <- NULL
bench_all <- NULL
for (nm in c("cycloess", "quantile")) {
  de <- process_proteomics(nadia_dia, norm_method = nm,
                           verbose = FALSE)$DEPs_results
  de$Species <- sp$PG.OrganismId[match(de$Protein.IDs, sp$PG.ProteinGroups)]
  opdea_all <- rbind(opdea_all, cbind(compute_opdea_metrics(de, ev), Assay = nm))
  bench_all <- rbind(bench_all, cbind(compute_benchmark_metrics(de, ev), Assay = nm))
}
ext <- bm_compute_extended_ranking(opdea_all, bench_all)
bm_plot_extended_ranking_bars(ext)
```
