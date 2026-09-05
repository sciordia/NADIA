# Compute extended ranking across 11 classification + OpDEA metrics

For each of the 11 metrics, averages across comparisons per Assay
(mean), then ranks methods in descending order (higher = better for all
metrics). Final rank = mean of the 11 individual ranks.

## Usage

``` r
bm_compute_extended_ranking(
  opdea_combined,
  bench_metrics_combined,
  metrics = .BM_EXTENDED_METRICS
)
```

## Arguments

- opdea_combined:

  data.frame from import_opdea_results()

- bench_metrics_combined:

  data.frame from import_benchmark_metrics()

- metrics:

  Character vector of metrics to include (default: all 11 extended
  metrics)

## Value

Named list with:

- extended_combined:

  Merged data.frame (long format)

- mean_aggregated:

  data.frame of mean metric values per Assay

- mean_ranking:

  data.frame with rank columns and rank_final

- n_methods:

  Number of unique methods

- n_comparisons:

  Number of unique comparisons

- metrics_used:

  Metrics actually used

## Details

Note: this extended ranking can drift away from the canonical OpDEA
ranking (`bm_compute_ranking`, based on nMCC/G_mean/pAUC). `Performance`
is collinear with F1 (so F1 effectively counts twice), and in spike-in
designs `Accuracy`/`NPV` are dominated by the background (TN \>\> ),
hence they have little discriminative power. Use it as a complementary
view, not as a substitute for the OpDEA ranking.

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
ext$mean_ranking[, c("Assay", "rank_final")]
#>      Assay rank_final
#> 1 cycloess       1.09
#> 2 quantile       1.91
```
