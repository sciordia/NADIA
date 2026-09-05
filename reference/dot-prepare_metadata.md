# Prepare metadata from proteomics_data object

Prepare metadata from proteomics_data object

## Usage

``` r
.prepare_metadata(preprocessing, covariate_df = NULL)
```

## Arguments

- preprocessing:

  A `proteomics_data` object, from any of the preprocess\_\*() functions

## Value

Data frame with columns: Column, Condition, Replicate
