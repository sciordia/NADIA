# Print a summary of a DIA-NN preprocessing result

Reports the number of runs, protein groups and conditions, so that the
shape of the experiment can be checked at a glance before processing it.

## Usage

``` r
# S3 method for class 'diann_data'
print(x, ...)
```

## Arguments

- x:

  A `diann_data` object, as returned by
  [`preprocess_diann()`](https://sciordia.github.io/NADIA/reference/preprocess_diann.md).

- ...:

  Ignored, present for compatibility with the `print` generic.

## Value

`x`, invisibly. Called for the summary it prints.

## Examples

``` r
diann <- preprocess_diann(
  system.file("extdata", "nadia_diann_report.tsv.gz", package = "NADIA"),
  condition_order = c("A", "B", "D"),
  verbose = FALSE)
print(diann)
#> Preprocessed DIA-NN data
#> ------------------------
#> Runs (metadata): 12 
#> Proteins (ID): 2000 
#> Proteins (QUANT): 2000 
#> 
#> Conditions: A, B, D 
```
