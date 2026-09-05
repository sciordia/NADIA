# Prepare PCA input with differential expression information

Prepare PCA input with differential expression information

## Usage

``` r
.prepare_pca_input(se, DEPs_results, assay_name, alpha = 0.05)
```

## Arguments

- se:

  SummarizedExperiment

- DEPs_results:

  Data frame with DE results

- assay_name:

  Assay name to use

- alpha:

  Significance threshold for sig_any

## Value

Data frame with intensities and adjP columns per comparison
