# Merge significance information from DEPs_results

Converts DEPs_results (long format) into wide format with one adjP\_\*
column per comparison, and computes sig_any.

## Usage

``` r
merge_significance_info(
  feature_ids,
  DEPs_results,
  assay_name,
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- feature_ids:

  Vector of feature IDs

- DEPs_results:

  DataFrame with columns Protein.IDs, adj.P.Val, Comparison, Assay

- assay_name:

  Name of the assay used to filter DEPs_results

- alpha:

  Significance threshold (default: 0.05)

- verbose:

  Print progress messages (default: TRUE)

## Value

DataFrame with FeatureID, adjP\_\*, sig_any
