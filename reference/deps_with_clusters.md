# Differential abundance with the Pattern Profiler cluster of each protein

Joins the differential-abundance results to the soft clustering, giving
the table a functional-enrichment or over-representation analysis needs:
the statistics that say which proteins moved, and the cluster that says
which pattern they follow.

## Usage

``` r
deps_with_clusters(
  DEPs_results,
  pattern_profiler,
  assignment = c("primary", "all"),
  min_membership = NULL,
  comparison = NULL,
  significant_only = FALSE
)
```

## Arguments

- DEPs_results:

  Differential-abundance results, the `DEPs_results` element of
  [`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md).
  Needs at least `Protein.IDs` and `Comparison`.

- pattern_profiler:

  Result of
  [`pattern_profiler_analysis()`](https://sciordia.github.io/NADIA/reference/pattern_profiler_analysis.md),
  or `NULL` for no clustering at all, in which case the three cluster
  columns are returned as `NA` and the shape of the table is unchanged.

- assignment:

  `"primary"` (default) for one row per protein and comparison, carrying
  the cluster of highest membership; `"all"` for one row per protein,
  comparison and cluster.

- min_membership:

  Optional extra membership threshold. It can only be stricter than the
  one the analysis ran with: pairs below that one were discarded when
  `long_output` was built and cannot be recovered here. A protein left
  with no cluster keeps its row, with `NA`.

- comparison:

  Optional character vector of comparisons to keep. The default, `NULL`,
  keeps them all.

- significant_only:

  Keep only the rows classified `Up` or `Down`. The classification is
  the one the analysis made, with its own `alpha` and `logFC_threshold`;
  nothing is recomputed here.

## Value

A data frame with the columns of `DEPs_results` followed by `Cluster`,
`Membership` and `ClusterRank` (1 is the protein's dominant cluster).

## Details

The join is a left join from the results, so **every protein that was
tested survives it**, with or without a cluster. That matters more than
it looks: the complete set of tested proteins is the background an
enrichment is computed against, and a protein with no cluster – one that
was never selected for the clustering, or that the clustering dropped –
belongs in that background rather than in the bin.

Clustering is fuzzy, so a protein can belong to several clusters at
once. `assignment` decides what to do about it. `"primary"` keeps only
the dominant pattern of each protein and returns exactly one row per
protein and comparison, which is what most enrichment tools expect.
`"all"` keeps every association above the threshold, which is the honest
representation of a soft clustering and the one to use when the analysis
downstream can weight by membership.

## See also

[`nadia_deps_with_clusters()`](https://sciordia.github.io/NADIA/reference/nadia_deps_with_clusters.md)
for the same table read from a `.nadia` file, and
[`pattern_profiler_analysis()`](https://sciordia.github.io/NADIA/reference/pattern_profiler_analysis.md)
for the clustering.

## Examples

``` r
data(nadia_dia)
res <- process_proteomics(nadia_dia, verbose = FALSE)

# Without a clustering the shape is the same and the clusters are NA.
fi <- deps_with_clusters(res$DEPs_results, NULL)
#> Warning: No Pattern Profiler result: the cluster columns are all NA.
head(fi)
#>         Protein.IDs Comparison Gene.Names       logFC      P.Value    adj.P.Val
#> 1        A0A024RBG1        B-A     NUDT4B  0.08293965 5.335637e-01 7.389228e-01
#> 2        A0A024RBG1        D-A     NUDT4B -0.01391145 9.160610e-01 9.511693e-01
#> 3        A0A024RBG1        D-B     NUDT4B -0.09685110 4.687369e-01 7.151013e-01
#> 4 A0A140T897;P02769        B-A        ALB  0.13439560 5.731060e-04 3.118509e-03
#> 5 A0A140T897;P02769        D-B        ALB  0.20481573 1.713786e-05 7.760613e-05
#> 6 A0A140T897;P02769        D-A        ALB  0.33921133 1.394707e-07 5.778484e-07
#>      Change         Assay MissGlobal MissComp MissCND1 MissCND2 Cluster
#> 1 No Change Impseqrob_min          0        0        0        0      NA
#> 2 No Change Impseqrob_min          0        0        0        0      NA
#> 3 No Change Impseqrob_min          0        0        0        0      NA
#> 4        Up Impseqrob_min          0        0        0        0      NA
#> 5        Up Impseqrob_min          0        0        0        0      NA
#> 6        Up Impseqrob_min          0        0        0        0      NA
#>   Membership ClusterRank
#> 1         NA          NA
#> 2         NA          NA
#> 3         NA          NA
#> 4         NA          NA
#> 5         NA          NA
#> 6         NA          NA

if (requireNamespace("Mfuzz", quietly = TRUE) &&
    requireNamespace("Biobase", quietly = TRUE) &&
    requireNamespace("e1071", quietly = TRUE)) {
    pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
                                    assay_name = "Impseqrob_min",
                                    verbose = FALSE)

    primary <- deps_with_clusters(res$DEPs_results, pp)
    nrow(primary) == nrow(res$DEPs_results)   # the background is intact

    soft <- deps_with_clusters(res$DEPs_results, pp, assignment = "all")
    table(soft$ClusterRank)
}
#> 
#>    1    2 
#> 2382  219 
```
