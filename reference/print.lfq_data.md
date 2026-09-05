# Print a summary of an LFQ preprocessing result

Reports the number of samples, protein groups and conditions, so that
the shape of the experiment can be checked at a glance before processing
it.

## Usage

``` r
# S3 method for class 'lfq_data'
print(x, ...)
```

## Arguments

- x:

  An `lfq_data` object, as returned by
  [`preprocess_lfq()`](https://sciordia.github.io/NADIA/reference/preprocess_lfq.md).

- ...:

  Ignored, present for compatibility with the `print` generic.

## Value

`x`, invisibly. Called for the summary it prints.

## Examples

``` r
lfq <- preprocess_lfq(
  system.file("extdata", "nadia_lfq_report.tsv.gz", package = "NADIA"),
  annot_path = system.file("extdata", "nadia_lfq_annotation.tsv",
                           package = "NADIA"),
  verbose = FALSE)
print(lfq)
#> Preprocessed LFQ (Proteome Discoverer) data
#> ------------------------------------------
#> Samples (metadata): 12 
#> Proteins (ID): 2000 
#> Proteins (QUANT): 2000 
#> 
#> Conditions: A, B, D 
```
