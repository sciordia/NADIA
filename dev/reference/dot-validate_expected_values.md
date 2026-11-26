# Validate expected_values

`expected_values` is the ground truth of the benchmark: it decides which
species are positives and which direction counts as a correct detection.
A malformed row does not produce an error further down, it produces a
plausible-looking benchmark, so every failure mode below is rejected
here.

## Usage

``` r
.validate_expected_values(ev)
```

## Arguments

- ev:

  Data frame with expected values

## Value

Invisible TRUE if valid, stops with error otherwise
