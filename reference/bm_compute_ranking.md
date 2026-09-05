# Compute OpDEA ranking across multiple methods

For each metric, aggregates across comparisons (mean and median), then
ranks methods in descending order (rank 1 = best). Final rank = average
of individual metric ranks.

## Usage

``` r
bm_compute_ranking(opdea_combined, metrics = .BM_METRICS)
```

## Arguments

- opdea_combined:

  data.frame with columns: Assay, Comparison, nMCC, G_mean, pAUC_001,
  pAUC_005, pAUC_010

- metrics:

  Character vector of metrics to include in ranking (default: all 5
  OpDEA metrics)

## Value

Named list with:

- mean_aggregated:

  data.frame of mean metric values per Assay

- median_aggregated:

  data.frame of median metric values per Assay

- mean_ranking:

  data.frame of ranks based on mean aggregation

- median_ranking:

  data.frame of ranks based on median aggregation

- opdea_combined:

  Input combined data.frame

- n_methods:

  Number of methods

- n_comparisons:

  Number of comparisons

- metrics_used:

  Metrics included in ranking

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
ranking$mean_ranking
#>      Assay rank_nMCC rank_G_mean rank_pAUC_001 rank_pAUC_005 rank_pAUC_010
#> 1 cycloess         1           1             1             1             1
#> 2 quantile         2           2             2             2             2
#>   rank_final
#> 1          1
#> 2          2
```
