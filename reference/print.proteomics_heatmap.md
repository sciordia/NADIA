# Draw a proteomics heatmap

Renders the object built by
[`proteomics_heatmap()`](https://sciordia.github.io/NADIA/reference/proteomics_heatmap.md)
on the current graphics device, adding the title stored in the object.
Printing is what actually draws the heatmap: building it does not.

## Usage

``` r
# S3 method for class 'proteomics_heatmap'
print(x, ...)
```

## Arguments

- x:

  A `proteomics_heatmap` object, as returned by
  [`proteomics_heatmap()`](https://sciordia.github.io/NADIA/reference/proteomics_heatmap.md).

- ...:

  Further arguments passed to
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).

## Value

The drawn `HeatmapList`, invisibly. Called for its side effect.

## Examples

``` r
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("tidyHeatmap", quietly = TRUE)) {
  data(nadia_dia)
  res <- process_proteomics(nadia_dia, verbose = FALSE)
  hm <- proteomics_heatmap(res$PCA_Input, mode = "any")

  # Drawing needs a graphics device; send it to a temporary file
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  print(hm)
  grDevices::dev.off()
  file.remove(f)
}
#> Registered S3 method overwritten by 'dendextend':
#>   method     from 
#>   rev.hclust vegan
#> tidyHeatmap says: (once per session) from release 1.7.0 the scaling is set to "none" by default. Please use scale = "row", "column" or "both" to apply scaling
#> tidyHeatmap says: If you use tidyHeatmap for scientific research, please cite: Mangiola, S. and Papenfuss, A.T., 2020. 'tidyHeatmap: an R package for modular heatmap production based on tidy principles.' Journal of Open Source Software. doi:10.21105/joss.02472.
#> This message is displayed once per session.
#> Warning: `when()` was deprecated in purrr 1.0.0.
#> ℹ Please use `if` instead.
#> ℹ The deprecated feature was likely used in the tidyHeatmap package.
#>   Please report the issue at
#>   <https://github.com/stemangiola/tidyHeatmap/issues>.
#> [1] TRUE
```
