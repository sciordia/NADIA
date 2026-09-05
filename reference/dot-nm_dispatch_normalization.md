# Dispatch a single normalization method on a log2 matrix

Internal helper that calls the appropriate `.norm_*()` function. For
Group A methods, converts log2 back to raw (2^x) before calling.

## Usage

``` r
.nm_dispatch_normalization(x_log2, method, method_args = list())
```

## Arguments

- x_log2:

  Numeric matrix in log2 scale (proteins x samples)

- method:

  Character scalar: normalization method name

- method_args:

  Named list of per-method arguments

## Value

Numeric matrix in log2 scale
