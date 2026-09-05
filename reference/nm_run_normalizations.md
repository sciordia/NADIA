# Run multiple normalization methods from a baseline assay

Takes a SummarizedExperiment with a log2-scale assay and applies each
requested normalization method, returning a new SE with one assay per
method.

## Usage

``` r
nm_run_normalizations(
  se,
  assay_name = "log2",
  methods = "all",
  method_args = list(),
  include_baseline = TRUE,
  verbose = TRUE
)
```

## Arguments

- se:

  SummarizedExperiment with at least one log2-scale assay.

- assay_name:

  Name of the baseline assay to normalize from. Default `"log2"`.

- methods:

  Character vector of method names, or `"all"` for all 12 benchmark
  methods. Default `"all"`.

- method_args:

  Named list of per-method arguments. E.g.
  `list(cycloess = list(method = "fast", span = 0.8))`.

- include_baseline:

  Logical. Include the baseline assay in the output SE. Default `TRUE`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

SummarizedExperiment with one assay per successfully normalized method
(plus baseline if `include_baseline = TRUE`).

## Examples

``` r
data(nadia_dia)
se <- nm_prepare_se(nadia_dia, verbose = FALSE)
se_bench <- nm_run_normalizations(se, methods = c("cycloess", "quantile", "MAD"),
                                  verbose = FALSE)
SummarizedExperiment::assayNames(se_bench)
#> [1] "log2"     "cycloess" "quantile" "MAD"     
```
