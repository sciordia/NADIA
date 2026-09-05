# Print a summary of a Spectronaut preprocessing result

Reports the number of runs, protein groups and conditions, so that the
shape of the experiment can be checked at a glance before processing it.

## Usage

``` r
# S3 method for class 'spectronaut_data'
print(x, ...)
```

## Arguments

- x:

  A `spectronaut_data` object, as returned by
  [`preprocess_spectronaut()`](https://sciordia.github.io/NADIA/reference/preprocess_spectronaut.md).

- ...:

  Ignored, present for compatibility with the `print` generic.

## Value

`x`, invisibly. Called for the summary it prints.

## Examples

``` r
data(nadia_dia)
print(nadia_dia)
#> Preprocessed Spectronaut data
#> -----------------------------
#> Runs (metadata): 12 
#> Proteins (ID): 2000 
#> Proteins (QUANT): 2000 
#> 
#> Conditions: A, B, D 
```
