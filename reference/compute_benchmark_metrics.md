# Compute benchmark metrics for all comparisons

Compute benchmark metrics for all comparisons

## Usage

``` r
compute_benchmark_metrics(
  de_res,
  ev,
  alpha = 0.05,
  lfc_thr = 0,
  p_col = "adj.P.Val",
  comparisons = NULL,
  assay = NULL,
  species_df = NULL
)
```

## Arguments

- de_res:

  Data frame with DE results

- ev:

  Data frame with expected values

- alpha:

  Significance threshold (default: 0.05)

- lfc_thr:

  Log fold-change threshold (default: 0)

- p_col:

  P-value column name (default: "adj.P.Val")

- comparisons:

  Comparisons to include (NULL = all)

- assay:

  Assay to filter (NULL = all)

- species_df:

  Data frame with Protein.IDs and Species columns (optional). If NULL,
  de_res must already contain a Species column.

## Value

Data frame with metrics per comparison

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
metrics[, c("Comparison", "TP", "FP", "Sensitivity", "Specificity", "F1")]
#>   Comparison  TP FP Sensitivity Specificity     F1
#> 1        B-A 469 27      0.6312      0.9785 0.7571
#> 2        D-A 693 60      0.9327      0.9522 0.9265
```
