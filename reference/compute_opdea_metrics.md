# Compute OpDEA metrics for all comparisons

Calculates additional benchmark metrics based on Peng et al. 2024
(Nature Comms, OpDEA framework):

- pAUC(0.01), pAUC(0.05), pAUC(0.1): Partial AUC with McClish correction

- nMCC: Normalized MCC = (MCC + 1) / 2

- G-mean: sqrt(Sensitivity \* Specificity)

## Usage

``` r
compute_opdea_metrics(
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

Data frame with columns: Comparison, nMCC, G_mean, pAUC_001, pAUC_005,
pAUC_010

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
compute_opdea_metrics(de, ev)
#>   Comparison   nMCC G_mean pAUC_001 pAUC_005 pAUC_010
#> 1        B-A 0.8411 0.7859   0.7677   0.8071   0.8211
#> 2        D-A 0.9412 0.9424   0.9125   0.9505   0.9614
```
