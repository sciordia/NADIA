# Compute OpDEA ranking separately for each comparison

For each comparison, ranks methods directly on their metric values (no
aggregation needed since there is one value per method per comparison).

## Usage

``` r
bm_compute_ranking_by_comparison(opdea_combined, metrics = .BM_METRICS)
```

## Arguments

- opdea_combined:

  data.frame with columns: Assay, Comparison, nMCC, G_mean, pAUC_001,
  pAUC_005, pAUC_010

- metrics:

  Character vector of metrics to include in ranking

## Value

Named list with:

- by_comparison:

  Named list where each key is a comparison and each value is a
  data.frame with metric values, ranks, and rank_final

- ranking_combined:

  data.frame with all comparisons in long format (includes Comparison
  column)

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
by_comp$by_comparison[["B-A"]]
#>      Assay   nMCC G_mean pAUC_001 pAUC_005 pAUC_010 rank_nMCC rank_G_mean
#> 1 quantile 0.8488 0.8040   0.7173   0.8046   0.8305         1           1
#> 2 cycloess 0.8411 0.7859   0.7677   0.8071   0.8211         2           2
#>   rank_pAUC_001 rank_pAUC_005 rank_pAUC_010 rank_final
#> 1             2             2             1        1.4
#> 2             1             1             2        1.6
```
