# NADIA 0.99.0

First version prepared for submission to Bioconductor. The project moves from a
collection of modules consumed with `source()` to an installable R package,
without any change to the results the pipeline produces.

## New

* Package structure: `DESCRIPTION`, `NAMESPACE` and help pages generated with
  roxygen2. The way to load the code is now `library(NADIA)`.
* 102 exported functions covering the whole pipeline: preprocessing of
  Spectronaut/DIA-NN and Proteome Discoverer reports (both TMT and label-free),
  normalization (13 methods), batch correction, imputation (19 methods),
  differential expression with limma or limpa, metrics and benchmarking, and
  interactive and static visualization.
* Example dataset `nadia_dia` and trimmed reports in `inst/extdata/`, with their
  provenance documented in `inst/scripts/`.

* Nine vignettes covering the pipeline end to end: `NADIA` (start here), plus
  `input-formats`, `missing-values`, `choosing-methods`, `benchmarking`,
  `batch-correction`, `pattern-profiler`, `visualization` and
  `results-and-export`. They replace `example_workflow.R`, which has been
  removed.
* A testthat suite of over 450 assertions, including one regression test per bug
  fixed in the 2026-07 code reviews.

## Changes

* **The missingness columns of `DEPs_results` are renamed, and there is a fourth
  one.** `MissingGlobal` was not global: it was the percentage of missing values
  over the replicates of the two conditions being compared, which coincides with
  the whole experiment only when there are exactly two conditions. Readers took
  the name at face value. The four columns are now, from the widest scope to the
  narrowest:

  | New | Old | Scope |
  |---|---|---|
  | `MissGlobal` | *(new)* | every sample in the experiment |
  | `MissComp` | `MissingGlobal` | the replicates of the two conditions compared |
  | `MissCND1` | `MissingPCT1` | the numerator condition |
  | `MissCND2` | `MissingPCT2` | the denominator condition |

  The old names are gone, with no aliases: code that reads `MissingGlobal` now
  fails loudly rather than returning the wrong figure quietly. The three renamed
  columns are identical value for value to the ones they replace — this is a
  rename, not a recomputation — and `MissGlobal` is new information that was not
  available anywhere before. The columns of the interactive tables are labelled
  `% Global`, `% Comp`, `% Cond 1` and `% Cond 2`, and each has its own numeric
  filter.
* Progress output goes through `message()` instead of `cat()`, so it can be
  silenced with `suppressMessages()` and redirected like any other condition.
  The four `print()` methods still use `cat()`, which is where it belongs. The
  text is unchanged except for thirteen strings that were still in Spanish.
* **Every random source is now under the caller's control.** `seed` arguments
  were added to `nm_compute_metrics()`, `normalization_metrics()` and
  `pattern_profiler_analysis()`. All three already had a seed internally, but no
  caller ever passed one, so the Hopkins statistic and — more importantly — the
  Mfuzz clustering could not be made reproducible on demand.
* **The PERMANOVA p-value is reproducible.** `vegan::adonis2()` derives it by
  permuting the group labels, so `PERMANOVA_pval`, a column that
  `normalization_metrics()` exports, changed between runs on the same data. It is
  now seeded from the same `seed` argument. `PERMANOVA_R2` was never affected.
* `quant_list_widget()` and `summary_list_widget()` no longer default `data` and
  `matrix_data` to files under `data-raw/`, which no user has. Both arguments are
  now required and the functions say so.
* `Depends: R (>= 4.6.0)`, and `biocViews` gains `BatchEffect`.
* Documentation notes that `missForest` takes no seed, and that the seed of `MLE`
  goes to the internal RNG of the `norm` package, which cannot be restored.

* **`max_na_prop` now defaults to `NULL`, which disables the pre-filter.** It
  used to default to 0.8, but the filter was only applied by `softHybrid` and by
  the single imputation methods; `combo`, the default, never applied it. The same
  argument therefore meant different things depending on `imp_method`. Passing a
  value still filters, now for every method alike. Runs that used the default
  `imp_method = "combo"` are unaffected; a run with a single method may now keep
  proteins that were previously discarded for having more than 80 % missing
  values — which in a package about missing values is the better default, since
  those are often the on/off cases.
* The "Export to Excel" button in the `*_widget()` tables no longer loads ExcelJS
  and PapaParse from a CDN. Both libraries (MIT) ship inside the package, so the
  button works offline and `saveWidget(selfcontained = TRUE)` embeds them. The
  exported `.xlsx` is byte-for-byte identical in content and formatting.
* Dependencies are declared in `Imports:` and `Suggests:`; the code no longer
  calls `library()`. Packages needed only for one specific method are checked at
  the point of use.
* Comments, documentation and user-facing messages are now in English.
* `process_proteomics()` and `pattern_profiler_analysis()` no longer write to
  disk unless given an output path: `export_dir` and `output_file` now default to
  `NULL`. `pattern_profiler_analysis()` returns `long_output` in its result list.
* Functions that call `set.seed()` restore the random number generator on exit,
  so they no longer affect the reproducibility of code run afterwards.
  `.nm_hopkins()` takes the seed as an argument.

## Bug fixes

* Eight duplicated definitions of colour and filtering helpers (`hex_to_rgba`,
  `darken_hex`, `normalize_hex`, `get_feature_ids` and others) spread across six
  modules are unified into one each. Because they all lived in the global
  environment, which copy was active depended on the order in which the modules
  had been loaded.
* `get_feature_ids()` now accepts `mode = "target"` and `mode = "specific"` as
  synonyms. Previously, loading `Heatmap_tidyHeatmap.R` after
  `PCA_Highcharts_Final.R` made `pca_highchart_list(modes = "specific")` abort.
* Filtering with `mode = "any"` can no longer return `NA` `FeatureID`s when the
  `sig_any` column contains missing values.
* `summary_list_widget()` works with LFQ and TMT metadata, which lack the
  identification-count columns specific to Spectronaut.
