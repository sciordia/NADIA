# NADIA 0.99.0

First version prepared for submission to Bioconductor. The project moves from a
collection of modules consumed with `source()` to an installable R package,
without any change to the results the pipeline produces.

## New

* Package structure: `DESCRIPTION`, `NAMESPACE` and help pages generated with
  roxygen2. The way to load the code is now `library(NADIA)`.
* 117 exported functions covering the whole pipeline: preprocessing of
  Spectronaut, DIA-NN and Proteome Discoverer reports (both TMT and label-free),
  normalization (13 methods), batch correction, imputation (20 methods),
  differential expression with limma or limpa, metrics and benchmarking, and
  interactive and static visualization.
* Example dataset `nadia_dia` and trimmed reports in `inst/extdata/`, with their
  provenance documented in `inst/scripts/`.
* Agent skill in `inst/skill/`, located with `nadia_skill_path()`: the pipeline
  written for a coding assistant, with the decision points and the guardrails
  the defaults invite.
* **`pattern_profiler_analysis()` no longer defaults `assay_name` to
  `"LoessCyc"`.** No pipeline built since the package conversion produces an
  assay with that name, so every call that omitted the argument failed. The
  default is now `NULL` and the assay is taken from `DEPs_results$Assay`, which
  records the one the differential abundance was computed on -- and which
  `assay_name` also filters, so the two have to agree anyway. When those results
  cover several assays, the function aborts and names them rather than guessing.
  Calls that already passed `assay_name` are unaffected.
* **`pattern_profiler_analysis(verbose = FALSE)` is now silent.** Its six
  helpers wrote their progress with bare `message()` calls, so a run asked to be
  quiet still reported the assay filter, the significance filter and the
  condition order. They take `verbose` now. `Mfuzz::filter.NA()` reports through
  `cat()`, which no argument of ours can reach, so its output is captured
  instead; and the Biobase startup banner, printed because NADIA attaches it for
  Mfuzz, is suppressed. The results are unchanged.
* **Batch correction is quiet too.** BERT reports through the `logging` package,
  which writes to stdout and honours no argument of NADIA's, so a corrected run
  printed fifteen `INFO::` lines whatever `verbose` said. They are captured now.
* Agent skill corrections, found by driving it through the TMT example from
  scratch: the `normalization_metrics()` call was documented without `methods`,
  which silently ranks the assays already in the object -- `raw`, `log2`, `BERT`,
  the imputed one -- instead of the twelve normalisation methods, and so would
  have justified a method choice with a comparison of pipeline stages; the join
  key of `covariate_df` (a `Column` column matching `metadata$Coding`, not
  `R.FileName`) was never stated; the import check rounded the missing
  percentage and reported 0 % for a dataset that had one missing cell; and
  nothing told the reader what to do when there is essentially nothing to
  impute.
* **Two imputation methods printed regardless of `verbose`.**
  `impute::impute.knn()` reports its recursive cluster split with `cat()` and
  `imputeLCMD::impute.MinProb()` prints its estimated sigma; neither is gated on
  anything, so both landed in quiet runs and in rendered vignettes. Captured in
  the wrappers, with a test that sweeps every method rather than the two known
  ones.
* A further skill correction: `imputation_metrics()` masks observed values, so
  its ranking measures MAR recovery and places every MNAR method last. Read as a
  straight ranking it says to drop the MNAR stage that the same dataset's
  missingness structure calls for. The skill now says to take only the MAR stage
  from that table, and to choose the MNAR stage from the intensity dependence or
  a spike-in.
* **`species_df` must map each protein exactly once.** The mapping is joined to
  the DE results on `Protein.IDs` with no uniqueness check, so a repeated
  identifier multiplied that protein's row in every comparison and inflated the
  counts every metric is built from -- one repeated protein took `classified_df`
  from 3994 to 3996 rows and TP from 469 to 470 on the example dataset, with no
  error and no warning. Exact repeats are now collapsed with a message, a
  protein claimed by two different species is an error (there is no way to know
  which side of the truth table it belongs to), and the join is asserted not to
  change the number of rows.

* **`expected_values` is validated as the ground truth it is.** A malformed
  benchmark specification did not produce an error further down, it produced a
  plausible-looking benchmark. `benchmarking_proteomics()` and the four other
  entry points now reject a non-numeric, `NA`, infinite or **zero**
  `expected_logFC` (a species expected not to change is background, and
  background is declared by being left out), blank or `NA` identifiers,
  duplicate `Comparison + Species` rows, and -- the quietest of them -- a
  **declared `Comparison` + `Species` pair with no proteins in the data**. A typo
  in an organism name used to leave those proteins in the background, so every
  correct detection of that species was counted as a false positive and the
  specificity collapsed with nothing in the output to say so. The check is per
  pair rather than global: a species present elsewhere in the results but absent
  from one of the comparisons that declares it leaves that comparison scored
  against a truth it does not contain, and if it was the only species declared
  there the comparison has no positives at all -- reporting `Sensitivity = NA`
  and `AUC = NA`, but `MCC = 0` and `nMCC = 0.5`, numbers rather than NAs, which
  then average silently into a multi-method ranking. A declared comparison that
  is absent altogether only warns: scoring a subset is legitimate.

* **Two diagnostic columns in `classified_df`**, also written to
  `benchmark_classified.tsv`. `is_significant` is the significance rule the
  classification actually applied, stored rather than recomputed so that no
  table or plot can drift from it. `direction_error` separates the two kinds of
  FN: a change missed altogether, and a change called significant with the sign
  inverted. Neither takes any part in a metric, and TP/FP/TN/FN are unchanged.

* **`preprocess_diann()`, the fourth input format.** DIA-NN protein-group
  matrices (`report.pg_matrix.tsv`) are wide — one row per protein group, one
  column per run — so they need their own reader rather than the Spectronaut
  one, which expects a long report. The result is the same `proteomics_data`
  object as the other three, and an example matrix ships in
  `inst/extdata/nadia_diann_report.tsv.gz`: the same twelve injections as
  `nadia_dia`, searched with DIA-NN instead of Spectronaut, which puts the
  missingness at 12.5 % against 8.1 % for the same runs.

  DIA-NN names its intensity columns after the raw file, which does not identify
  the experimental design, so the design is declared rather than guessed: either
  rename the columns to the `Abundance: <Condition>_<Replicate>` convention that
  `preprocess_tmt()` already reads, or pass a sample sheet through `annot_path`
  and leave the report untouched. `R.FileName` keeps the original header either
  way.
  A protein-group matrix carries no per-sample metrics at all; what DIA-NN does
  report is mapped onto the wide contract, and what it does not is left out
  rather than filled with `NA`, except for the two families the interactive
  protein and quantification tables need in order to detect their sample
  columns.

* **`imp_method = "halfmin"`**, the half-minimum imputation that DIA-NN uses for
  its own differential analysis: "replaces each missing value with half the
  observed minimum across all runs". It can be used on its own or as the MNAR
  stage of `combo` and `softHybrid` (`mnar_method = "halfmin"`), and it takes part
  in the `imputation_metrics()` benchmark. Since the assay is on the log2 scale,
  halving the intensity means **one unit below** the global observed minimum, not
  half of the log2 value — `log2(m / 2)` is `log2(m) - 1`. `combo` still defaults
  to `mnar_method = "min"`; nothing existing changes.

* **The `.nadia` format: one file for a whole analysis.** `write_nadia()` writes
  the preprocessing tables, the processing results and, optionally, the Pattern
  Profiler output to a single [DuckDB](https://duckdb.org) database, together
  with the parameters and the provenance of the run. `read_nadia()` reads it
  back, and `nadia_result()` and `nadia_preprocessing()` return objects the rest
  of the package plots directly — the round trip is `identical()`, column order,
  integer columns and ordered factors included.

  Nothing is stored twice. The `log2` assay is exactly `log2` of the raw one and
  is not stored; `PCA_Input` is the imputed assay joined to the results and is
  not stored either; and of the imputed assay only the cells actually filled are
  kept, each tagged with the branch that filled it — so `imputed_values` is both
  the imputed data and the MAR/MNAR mask. On the example dataset the file is
  4.1 MB against 6.8 MB of TSV and 2.7 MB of Parquet, in one file instead of
  eleven.

  The reconstruction lives in the file as SQL views rather than in R code, so a
  single definition serves every client; with backward-compatible storage, that
  includes DuckDB-WASM in a browser, with no R installed. `nadia_connect()`
  opens the file for SQL, `nadia_tables()` lists it without loading it, and
  `nadia_export_parquet()` is the way out to an open format. Needs `duckdb` and
  `DBI`, both in `Suggests`.
* **`deps_with_clusters()` and `nadia_deps_with_clusters()`: the input a
  functional analysis needs.** Differential abundance says which proteins moved
  and Pattern Profiler says which pattern they follow; an enrichment needs both
  in one table, and joining them by hand is easy to get wrong — `DEPs_results`
  is long by comparison, `long_output` is long by cluster, and merging them
  without saying what to do about that produces a silent cartesian product.
  The two functions return the columns of `DEPs_results` followed by `Cluster`,
  `Membership` and `ClusterRank`, either from the two results in memory or from
  a `.nadia` file, where the join is the view `v_deps_pattern_profiler`. They
  return identical tables, which the tests assert rather than assume.

  The join is a left join from the results and stays one: every protein that
  was tested survives it, with or without a cluster. That set is the background
  an enrichment is measured against, and quietly dropping part of it changes
  every result computed from it. `assignment = "primary"` gives one row per
  protein and comparison, carrying the dominant cluster; `assignment = "all"`
  keeps every membership above the threshold, which is what a soft clustering
  actually says.

  Feeding that needed one addition to the file. `min_membership` is applied when
  `long_output` is built, so a protein whose highest membership falls below it
  has no row at all — and the highest membership of a protein is only
  guaranteed to be at least `1 / optimal_c`, so with enough clusters this
  happens. Stored that way, a protein that *was* clustered is indistinguishable
  from one that never entered the clustering. The new `pp_assignment` table
  holds the unfiltered hard assignment, one row per clustered protein, and
  `"primary"` now returns a cluster for every one of them. On the example
  dataset it costs 794 rows.

  The addition is backward compatible in both directions: a file written
  earlier has neither the table nor the view and reads without them, and a file
  written now is read by an earlier NADIA that ignores them. The schema version
  is unchanged, deliberately — adding an object is not changing the shape of
  one.
* Nine vignettes covering the pipeline end to end: `NADIA` (start here), plus
  `input-formats`, `missing-values`, `choosing-methods`, `benchmarking`,
  `batch-correction`, `pattern-profiler`, `visualization` and
  `results-and-export`. They replace `example_workflow.R`, which has been
  removed.
* A testthat suite of over 450 assertions, including one regression test per bug
  fixed in the 2026-07 code reviews.

## Changes

* **The label-free example is now the same experiment as the DIA ones.** It used
  to be an unrelated label-free run of six samples in two conditions (WT and
  MUT). `inst/extdata/nadia_lfq_report.tsv.gz` and its annotation now hold the
  same twelve injections as `nadia_dia` and `nadia_diann` — conditions A, B and D
  with four replicates each — acquired label-free and searched with Proteome
  Discoverer 3.3. Comparing the three example formats is now comparing
  acquisition and search methods rather than comparing experiments: 8.0 % missing
  values against 8.1 % for Spectronaut and 12.5 % for DIA-NN on the same runs.
  Anything reading the example gets three contrasts (`B-A`, `D-A`, `D-B`) where
  it used to get one (`MUT-WT`).

* **The three wide readers declare their design the same way.** Proteome
  Discoverer LFQ used to be the odd one out: it *required* a separate annotation
  file, while TMT and DIA-NN read the `Abundance: <Condition>_<Replicate>`
  column suffixes. Those suffixes were there all along -- the sheet the package
  ships is a table restating them -- so `preprocess_lfq()` now reads them too,
  and `annot_path` is optional in all three.

  `sample_names` is retired from `preprocess_diann()` in favour of the same
  `annot_path`. A vector of labels has to be positional, and nothing guarantees
  the column order of an export: Proteome Discoverer does not even write the
  four LFQ column families in the same order as each other, so a re-export that
  reordered runs would have relabelled them all without any way to notice. A
  sheet is matched by name, so it either matches or it aborts.

  The sheet needs `Column` and `Condition`, where `Column` is the sample
  identifier as it appears in the report. Replicates come from a trailing
  `_<digits>` when the identifier has one, so a sheet restating the suffixes
  gives a result `identical()` to supplying no sheet at all. Anything else in
  the sheet is ignored, including the `Experiment` column the LFQ annotation
  used to carry: batch and experiment structure belongs in `covariate_df` of
  `process_proteomics()`. The `R.Experiment` column is gone from `metadata`,
  where nothing ever read it.

* **The Spectronaut export schema ships with the package.** Reading a Spectronaut
  report requires it to have been exported with the right columns, which until
  now meant assembling that list by hand from the documentation.
  `system.file("extdata", "NADIA_Report.rs", package = "NADIA")` is the schema
  itself: import it into Spectronaut and export with it, and the report has the
  shape `preprocess_spectronaut()` expects.

* **The TMT example drops condition C, so all four example reports share one
  design.** The TMTpro report is a separate experiment from the other three, but
  it reproduces the same three-proteome spike-in — *E. coli* rising from A to D
  and yeast falling, with human and bovine flat — so its condition labels already
  meant what they mean elsewhere. It now ships conditions A, B and D of eight
  replicates (2,000 × 43), and `condition_order = c("A", "B", "D")` replaces
  `c("A", "B", "C", "D")` in the examples and the vignette.

  The two TMT mixes are deliberately kept. They are the only real batch structure
  among the examples and a balanced one — replicates 1–4 are the first mix and 5–8
  the second, bridged by the four `IS` channels — with the mix accounting for 67 %
  of the variance on PC1 and the condition for 29 % on PC2, orthogonally. A
  Proteome Discoverer report carries no mix column, so that mapping is now stated
  in the vignette and in the help page rather than left implicit.

* **`preprocess_lfq()` and `preprocess_tmt()` read the Proteome Discoverer 3.3
  column names.** Version 3.3 renamed `Exp. q-value: Combined` to
  `Exp. Protein q-value: Combined`. Both readers resolve a column by trying a
  list of known spellings and fall back to `NA` when none matches, so a report
  from the new version used to produce an `Exp.q.value` column that was entirely
  `NA`, silently. The new spelling is now among the candidates in both readers;
  older exports keep working unchanged.

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

* **A Pattern Profiler added to a preprocessing-only `.nadia` file was written
  but never readable.** `.db_create_views()` returned early when the file held
  no processing step, and the block that builds `v_pattern_profiler` sat after
  that return. So `nadia_add_pattern_profiler()` on a file written with
  `write_nadia(f, preprocessing)` alone wrote `pp_membership` and `pp_profile`,
  set `has_pattern_profiler` to `TRUE`, and created no view: the metadata said
  the clustering was there and `nadia_pattern_profiler()` returned `NULL`, with
  nothing to say why. The view is now built before the gate.

* **The "Significant Proteins" plot left out the proteins it was best placed to
  expose.** `benchmark_signif_bars_gg()` filtered on `predicted == 1`, which is
  significance *and* the expected direction. A spike-in found significant with
  the sign inverted is an FN and so carries `predicted = 0`, so it disappeared
  from a chart whose title promises every significant protein. On the example
  dataset nine yeast proteins were missing, and they were missing from the UP
  facet -- the plot showed no yeast protein rising in any comparison when nine
  did. The filter now reads the new `is_significant` column.

* **The significant-protein summary ignored `lfc_thr`.**
  `.summarize_significant_proteins()` counted on `p <= alpha` alone, while the
  classification behind every metric also required `abs(logFC) >= lfc_thr`, so
  the table and the metrics in the same output disagreed as soon as a
  fold-change threshold was in force: 1922 significant proteins reported against
  1857 classified on the example dataset at `lfc_thr = 0.3`. It now reads
  `is_significant` instead of recomputing the rule, which makes the two
  structurally unable to diverge. With `lfc_thr = 0`, the default, the table is
  unchanged.

* **`batch_covariates` now works with categorical metadata.** Protecting a
  biological variable during batch correction aborted for every non-numeric
  column, which is the usual case: a condition is stored as text. BERT performs
  no encoding of its own -- it collects the `Cov_*` columns and passes them
  straight to `sva::ComBat(mod = )` and `limma::removeBatchEffect(design = )` --
  and NADIA handed it the column unchanged. ComBat coerced the text to double,
  produced `NA`s and failed with `NA/NaN/Inf in foreign function call`, while
  limma rejected it with `design must be a numeric matrix`. Categorical
  covariates are now encoded as treatment-contrast indicators before the call.
  Numeric covariates are still passed through untouched, so a continuous
  covariate keeps its meaning and results that already worked are unchanged.

* `process_proteomics()` no longer emits a `max()` warning for every protein
  quantified in no sample by the search engine. The unique-peptide count is the
  row maximum of the per-sample peptide columns, and a row of all `NA` made
  `max(na.rm = TRUE)` warn and return `-Inf`. The value was already replaced by
  zero on the next line, so the warning was noise; a label-free report produces
  one per protein detected by the feature detector but never scored.

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
