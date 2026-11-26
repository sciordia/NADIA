# Generate comparisons for differential analysis

Generate comparisons for differential analysis

## Usage

``` r
.specify_comparisons(se, condition_column = NULL, control = NULL)
```

## Arguments

- se:

  SummarizedExperiment

- condition_column:

  Condition column (if NULL, uses SE metadata)

- control:

  Control condition. If NULL, generates all pairwise comparisons

## Value

Factor with comparisons in "Treatment-Control" format
