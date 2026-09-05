# Build the Colour Palette for the Heatmap Values

Build the Colour Palette for the Heatmap Values

## Usage

``` r
get_heatmap_palette(palette = NULL, n = 11, reverse = FALSE)
```

## Arguments

- palette:

  Palette specification:

  - NULL: use the default palette (RdBu)

  - Vector of colours: use those colours directly

  - String "paletteer::" (e.g. "viridis::viridis"): use paletteer

  - String "brewer:" (e.g. "brewer:RdBu"): use RColorBrewer

- n:

  Number of colours to generate

- reverse:

  Reverse the palette (default: FALSE)

## Value

Vector of colours, or a colorRamp2 function
