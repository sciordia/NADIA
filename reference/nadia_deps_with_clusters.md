# Differential abundance with clusters, from a `.nadia` file

The file's own answer to
[`deps_with_clusters()`](https://sciordia.github.io/NADIA/reference/deps_with_clusters.md).
The join lives in the file as the view `v_deps_pattern_profiler`, so
this reads a table rather than building one, and a client that is not R
gets the same table from the same SQL. The filters are the ones
[`deps_with_clusters()`](https://sciordia.github.io/NADIA/reference/deps_with_clusters.md)
applies, and the two functions return identical results for the same
arguments.

## Usage

``` r
nadia_deps_with_clusters(
  db,
  assignment = c("primary", "all"),
  min_membership = NULL,
  comparison = NULL,
  significant_only = FALSE
)
```

## Arguments

- db:

  A `nadia_db` from
  [`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
  or the path to a `.nadia` file.

- assignment:

  `"primary"` (default) for one row per protein and comparison, carrying
  the cluster of highest membership; `"all"` for one row per protein,
  comparison and cluster.

- min_membership:

  Optional extra membership threshold. It can only be stricter than the
  one the analysis ran with.

- comparison:

  Optional character vector of comparisons to keep.

- significant_only:

  Keep only the rows classified `Up` or `Down`.

## Value

A data frame with the columns of `DEPs_results` followed by `Cluster`,
`Membership` and `ClusterRank`, or `NULL` if the file holds no
differential-abundance results. A file with results but no clustering
returns the same shape with the three cluster columns `NA`.

## Details

The table is read on demand rather than by
[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
because it is as long as the number of proteins times comparisons times
clusters and there is no reason to carry it around unasked.

## See also

[`deps_with_clusters()`](https://sciordia.github.io/NADIA/reference/deps_with_clusters.md),
[`nadia_add_pattern_profiler()`](https://sciordia.github.io/NADIA/reference/nadia_add_pattern_profiler.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    f <- file.path(tempdir(), "functional.nadia")
    write_nadia(f, nadia_dia, res, verbose = FALSE)

    head(nadia_deps_with_clusters(f))

    unlink(f)
}
#> Warning: 'functional.nadia' holds no Pattern Profiler view: the cluster columns are all NA. Add a clustering with nadia_add_pattern_profiler().
```
