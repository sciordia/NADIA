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

## Changes

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
