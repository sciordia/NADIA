# halfmin: half of the global observed minimum (the DIA-NN recipe)

DIA-NN describes it as replacing each missing value with half the
observed minimum across all runs. The matrix here is on the log2 scale,
so halving the intensity means subtracting 1, not dividing by 2: log2(m
/ 2) == log2(m) - 1. Dividing the log2 value would give the square root
of the minimum intensity, eight times lower than intended on a typical
dataset, and its magnitude would depend on where zero happens to fall on
the log scale.

## Usage

``` r
.imp_halfmin(x, args = list())
```
