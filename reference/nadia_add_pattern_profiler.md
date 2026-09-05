# Add a Pattern Profiler run to an existing `.nadia` file

The clustering is an optional extra step, so it is written separately:
this adds it to a file that already holds the analysis, and refuses to
create one from scratch.

## Usage

``` r
nadia_add_pattern_profiler(
  file,
  pattern_profiler,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- file:

  Path to an existing `.nadia` file.

- pattern_profiler:

  Object from
  [`pattern_profiler_analysis()`](https://sciordia.github.io/NADIA/reference/pattern_profiler_analysis.md).

- overwrite:

  Replace a clustering already in the file.

- verbose:

  Print progress.

## Value

The path, invisibly.

## See also

[`write_nadia()`](https://sciordia.github.io/NADIA/reference/write_nadia.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    res <- process_proteomics(nadia_dia, verbose = FALSE)
    f <- file.path(tempdir(), "pp.nadia")
    write_nadia(f, nadia_dia, res, verbose = FALSE)

    # pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results)
    # nadia_add_pattern_profiler(f, pp)

    unlink(f)
}
```
