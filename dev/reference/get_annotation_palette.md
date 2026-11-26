# Build the Colour Palette for Categorical Annotations

Build the Colour Palette for Categorical Annotations

## Usage

``` r
get_annotation_palette(levels, palette = NULL)
```

## Arguments

- levels:

  Vector of levels (categories)

- palette:

  Colour palette: `NULL` for the default; `"pkg::name"` for a palette
  from paletteer; `"brewer:Name"`, or the bare name, for an RColorBrewer
  palette; or a vector of colours, unnamed or named after the groups. A
  palette with fewer colours than levels is interpolated.

## Value

Named vector of colours
