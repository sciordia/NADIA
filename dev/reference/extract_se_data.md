# Extract data from a SummarizedExperiment

Retrieves the intensity matrix, the feature IDs and the sample metadata
from a SummarizedExperiment object.

## Usage

``` r
extract_se_data(se_proc, assay_name = NULL, verbose = TRUE)
```

## Arguments

- se_proc:

  SummarizedExperiment object

- assay_name:

  Name of the assay to use (NULL = first available)

- verbose:

  Print progress messages (default: TRUE)

## Value

List with:

- intensity_matrix: features x samples matrix

- feature_ids: vector of protein IDs

- sample_metadata: DataFrame with SampleID, Condition, etc.
