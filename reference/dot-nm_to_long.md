# Reshape SE assays to long data.frame for ggplot2

Reshape SE assays to long data.frame for ggplot2

## Usage

``` r
.nm_to_long(se, assay_names = NULL, condition_col = "Condition")
```

## Arguments

- se:

  SummarizedExperiment object

- assay_names:

  Character vector of assay names to include. NULL = all.

- condition_col:

  Column name in colData with condition labels.

## Value

data.frame with columns: Method, Sample, Condition, Protein, Value
