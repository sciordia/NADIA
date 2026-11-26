# Detect elbow point via maximum perpendicular distance to LOESS fit

Given missing_rate (x) and mean_intensity (y), fits a LOESS curve and
finds the point of maximum perpendicular distance to the line connecting
the first and last fitted points.

## Usage

``` r
.detect_elbow(missing_rate, mean_intensity)
```

## Arguments

- missing_rate:

  Numeric vector of per-protein missing rates

- mean_intensity:

  Numeric vector of per-protein mean intensities (non-NA)

## Value

List with r0 (elbow missing rate) and x0 (elbow intensity)
