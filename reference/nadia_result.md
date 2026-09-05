# Rebuild the processing result from a `.nadia` file

Returns a `proteomics_result` with everything
[`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md)
returns, the `SummarizedExperiment` included, so that every plotting
function in the package works from the file directly:

## Usage

``` r
nadia_result(db)
```

## Arguments

- db:

  A `nadia_db` from
  [`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md).

## Value

A `proteomics_result` object, or an error if the file holds only the
preprocessing step.

## Details

    db  <- read_nadia("experiment.nadia")
    res <- nadia_result(db)
    volcano_highchart_list(res$DEPs_results)
    pca_highchart_list(res$PCA_Input)

## See also

[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
[`nadia_preprocessing()`](https://sciordia.github.io/NADIA/reference/nadia_preprocessing.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    res <- process_proteomics(nadia_dia, verbose = FALSE)
    f <- file.path(tempdir(), "result.nadia")
    write_nadia(f, nadia_dia, res, verbose = FALSE)

    back <- nadia_result(read_nadia(f))
    identical(back$DEPs_results, res$DEPs_results)
    class(volcano_highchart_list(back$DEPs_results)[[1]])

    unlink(f)
}
```
