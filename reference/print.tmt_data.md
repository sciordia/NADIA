# Print a summary of a TMT preprocessing result

Reports the number of channels, protein groups and conditions, so that
the shape of the experiment can be checked at a glance before processing
it.

## Usage

``` r
# S3 method for class 'tmt_data'
print(x, ...)
```

## Arguments

- x:

  A `tmt_data` object, as returned by
  [`preprocess_tmt()`](https://sciordia.github.io/NADIA/reference/preprocess_tmt.md).

- ...:

  Ignored, present for compatibility with the `print` generic.

## Value

`x`, invisibly. Called for the summary it prints.

## Examples

``` r
tmt <- preprocess_tmt(
  system.file("extdata", "nadia_tmt_report.tsv.gz", package = "NADIA"),
  condition_order = c("A", "B", "D"),
  verbose = FALSE)
print(tmt)
#> Preprocessed TMT (Proteome Discoverer) data
#> ------------------------------------------
#> Channels (metadata): 24 
#> Proteins (ID): 2000 
#> Proteins (QUANT): 2000 
#> 
#> Conditions: A, B, D 
```
