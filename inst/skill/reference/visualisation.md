# Visualisation: which figure answers which question

Every function takes tables or data frames, so they work with the output of
`process_proteomics()` or with results produced elsewhere.

| Question | Function | Input |
|---|---|---|
| Which proteins change, and by how much? | `volcano_highchart_list()` | `DEPs_results` |
| Did normalisation work? | `boxplot_highchart_list()` | `BoxPlot_Input` |
| Do samples group by condition? | `pca_highchart_list()` | `PCA_Input` |
| What is the pattern across conditions? | `proteomics_heatmap_list()` | `PCA_Input` |
| Which proteins behave alike? | `pattern_profiler_analysis()` + the cluster charts | `se_proc` + `DEPs_results` |
| What are the numbers, exactly? | the `*_reactable()` / `*_widget()` tables | see `reference/export.md` |

Highcharts functions return **named lists**, one entry per comparison or mode.
Pick one to display or save:

```r
volcanoes <- volcano_highchart_list(res$DEPs_results, comparisons = "B-A",
                                    alpha = 0.05, lfc_thr = 1)
volcanoes[["B-A"]]
```

## The interactive ones

```r
# Volcano
volcano_highchart_list(de_res, ain, comparisons, lfc_thr, alpha, p_col,
                       highlight_genes, show_top_genes, palette, title)

# Boxplots: one box per sample, faceted by assay -> shows what normalisation did
boxplot_highchart_list(data, assays, color_by, group_order, palette, horizontal)

# PCA, with confidence ellipses or hulls
pca_highchart_list(pca_input, modes, alpha, color_by, addEllipses,
                   ellipse_type, ellipse_level, filter_samples_to_comparison)
```

`boxplot_highchart_list()` is the honest check on normalisation: pass both the
`log2` and the normalised assay and look at whether the distributions align.
If they do not, the method did not do its job.

`pca_highchart_list(modes = "specific")` restricts each plot to the samples of
its comparison, which is what you want when conditions differ in n.

## The static ones

```r
proteomics_heatmap_list(
  data,                       # PCA_Input
  modes,                      # e.g. "significant"
  alpha = 0.05,
  scale_data = TRUE,          # row z-score
  split_by_condition = TRUE,
  export_path = "figures", export_modes = "png"
)
```

Needs `tidyHeatmap`/`ComplexHeatmap` (`Suggests`). ggplot2 figures also come
out of the metrics and benchmarking orchestrators.

## Pattern Profiler: proteins with similar profiles

Soft (fuzzy) clustering with Mfuzz. **A separate step, not part of
`process_proteomics()`.**

```r
pp <- pattern_profiler_analysis(
  res$se_proc,
  res$DEPs_results,
  # assay_name  = NULL,            # resolved from DEPs_results$Assay
  filter_mode = "any",             # "any" | "all" | "specific"
  alpha       = 0.05,
  c_range     = 2:10,
  auto_select_c   = TRUE,
  selection_method = "xb",         # "xb" | "consensus" | "elbow"
  min_membership  = 0.25,
  seed = 123,
  verbose = FALSE
)

pp$optimal_c
cluster_profile_highchart_list(pp$long_output)
cluster_centroids_highchart(pp$long_output)
```

**Three things to know:**

1. `assay_name` selects the matrix to cluster **and** filters `DEPs_results` by
   its `Assay` column, so the two must agree. Left at `NULL` it is taken from
   `DEPs_results$Assay`. Pass it to cluster a different assay, or when the
   results cover several — NADIA then aborts rather than guessing.
2. `filter_mode = "all"` means **all features**, i.e. no significance filter at
   all. It does *not* mean "significant in all comparisons". `"any"` is the one
   that keeps proteins significant in at least one comparison.
3. Mfuzz needs `Biobase` and `e1071` **attached**. The function does this
   itself and releases them on exit; a hand-written Mfuzz call will not.

Membership is soft: a protein belongs to every cluster with a weight.
`min_membership` filters what reaches `long_output`, so a protein whose best
membership falls below it disappears from the long table entirely.

## Joining clusters to the statistics

For enrichment or functional analysis:

```r
fi <- deps_with_clusters(res$DEPs_results, pp, assignment = "primary")
```

- `assignment = "primary"` — one row per protein and comparison, the dominant
  cluster (`ClusterRank == 1`).
- `assignment = "all"` — every retained soft association; more rows.

It is a **left join from `DEPs_results` and must stay one**: proteins with no
cluster keep `NA` and stay in the table, because the full set of tested
proteins is the enrichment background. Do not drop them.

## Colour palettes

Built-in names, RColorBrewer and paletteer sources are all accepted through
`palette`. When the number of clusters exceeds the 15-colour default palette,
NADIA interpolates the extra colours rather than failing, so a large `c` needs
no special handling.

## A note on Highcharts licensing

Interactive figures use `highcharter`, an interface to the Highcharts
JavaScript library, which has its own licence terms. NADIA's GPL licence grants
no right to use Highcharts. For commercial or institutional use, check
<https://www.highcharts.com/license>. The whole statistical pipeline runs
without producing a single Highcharts figure — ggplot2 and ComplexHeatmap
cover the static output.
