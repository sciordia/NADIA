# Build the Colour Palette for the Heatmap Values

Build the Colour Palette for the Heatmap Values

## Usage

``` r
get_heatmap_palette(palette = NULL, n = 11, reverse = FALSE)
```

## Arguments

- palette:

  Palette for the values, used as a colour ramp: `NULL` for the default
  (blue-white-red); `"pkg::name"` for a palette from paletteer;
  `"brewer:Name"`, or the bare name, for an RColorBrewer palette (up to
  `n` colours); or a vector of colours.

- n:

  Number of colours to generate

- reverse:

  Reverse the palette (default: FALSE)

## Value

Vector of colours, or a colorRamp2 function
