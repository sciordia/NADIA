# Build the clustering matrix aggregated by condition

Build the clustering matrix aggregated by condition

## Usage

``` r
build_clustering_matrix(
  intensity_matrix,
  sample_metadata,
  selected_features,
  condition_order = NULL,
  aggregate = c("median", "mean"),
  verbose = TRUE
)
```

## Arguments

- intensity_matrix:

  Intensity matrix (features x samples)

- sample_metadata:

  DataFrame with SampleID, Condition

- selected_features:

  Vector of features to include

- condition_order:

  Condition order (NULL = alphabetical order)

- aggregate:

  Aggregation method: "median" or "mean"

- verbose:

  Print progress messages (default: TRUE)

## Value

Matrix of aggregated intensities (features x conditions)
