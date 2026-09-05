# Convert zero values to NA

Replaces 0 values in a matrix or data frame with NA. Useful for DIA-NN
data where 0 indicates non-detection.

## Usage

``` r
.zero_to_missing(data)
```

## Arguments

- data:

  Numeric matrix or data frame

## Value

Object of same type with 0 converted to NA
