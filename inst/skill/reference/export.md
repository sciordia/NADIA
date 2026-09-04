# Export and archiving

Four routes, for different purposes.

| Route | Best for | Main limitation |
|---|---|---|
| TSV | inspection in a text editor or spreadsheet | several files to keep together; types inferred on read |
| Parquet | typed exchange of large tables | needs Parquet-aware software |
| `.nadia` | preserving and restoring a complete analysis | needs NADIA or a DuckDB client |
| Interactive HTML | exploring or presenting without R | may need a dependency folder |

## TSV and Parquet, straight from the pipeline

Export happens only when `export_dir` is supplied. With the default `NULL`,
nothing is written.

```r
process_proteomics(
  prep,
  export_dir    = "results/tables",
  export_format = "tsv",        # "tsv" | "parquet" | "both"

  export_normalized = TRUE,     # the normalised matrix
  export_imputed    = TRUE,     # the final imputed matrix
  export_volcano    = TRUE,     # the DE results
  export_boxplot    = TRUE,     # long-format intensities
  export_pca        = TRUE      # long-format + significance flags
)
```

Five datasets by default. File names carry the normalisation and imputation
methods, so runs do not get mixed up:

```
matrix_log2_cycloess.tsv
matrix_log2_cycloess_Impseqrob_min.tsv
VolcanoPlot_Input_cycloess_Impseqrob_min.tsv
BoxPlot_Input_cycloess_Impseqrob_min.tsv
PCA_Input_cycloess_Impseqrob_min.tsv
```

Parquet through this route needs `arrow`.

## The `.nadia` archive

One DuckDB file holding the whole analysis: preprocessing, processing,
parameters, provenance, and optionally the clustering.

```r
write_nadia("analysis.nadia",
            preprocessing    = prep,
            result           = res,
            pattern_profiler = pp,     # optional
            overwrite        = FALSE,
            verbose          = FALSE)

db  <- read_nadia("analysis.nadia")
res2 <- nadia_result(db)            # reconstructs the process_proteomics() output
prep2 <- nadia_preprocessing(db)    # reconstructs the proteomics_data object
```

The round trip is `identical()`, not merely `all.equal()`: the file records
column order, R classes and factor levels.

Nothing derivable is stored twice. The `log2` assay and `PCA_Input` are
reconstructed by SQL views; for imputation methods that leave observed values
untouched, only the filled cells are kept — which makes that table double as
the MAR/MNAR mask.

### Inspecting and querying

```r
nadia_tables("analysis.nadia")      # tables and views, without loading

con <- nadia_connect("analysis.nadia")
DBI::dbGetQuery(con,
  'SELECT "Comparison", "Change", COUNT(*) AS n
     FROM v_deps_results GROUP BY 1, 2 ORDER BY 1, 2')
DBI::dbDisconnect(con, shutdown = TRUE)
```

Canonical views: `v_metadata`, `v_protein_id`, `v_protein_quant`,
`v_deps_results`, `v_boxplot_input`, `v_pca_input`, `v_matrix_norm`,
`v_matrix_imputed`, and — when clustering is stored — `v_pattern_profiler` and
`v_deps_pattern_profiler`.

Adding clustering afterwards, and reading the joined table:

```r
nadia_add_pattern_profiler("analysis.nadia", pattern_profiler = pp)
nadia_deps_with_clusters(read_nadia("analysis.nadia"), assignment = "primary")
```

### Extracting to Parquet

```r
nadia_export_parquet(db, "results/parquet")
```

DuckDB writes the files itself, so `arrow` is **not** needed on this route. It
exports the canonical views plus the metadata tables (`nadia_meta`,
`nadia_calls`, `nadia_packages`, `nadia_source_files`, `parameters`).

### One hard rule

**Never keep an open `.nadia` file inside iCloud Drive, Dropbox or OneDrive.**
The sync client will corrupt an open DuckDB database. Write and use it in a
local directory, close it, and copy the finished file across.

`duckdb` and `DBI` are in `Suggests`.

## Interactive HTML

Four interactive tables. `*_reactable()` returns the table widget alone;
`*_widget()` returns a complete page with a toolbar and Excel export. Only
results and protein identification have both variants.

| Content | Page | Table only |
|---|---|---|
| Sample metadata | `summary_list_widget()` | — |
| Protein annotations | `protein_list_widget()` | `protein_list_reactable()` |
| Quantification values | `quant_list_widget()` | — |
| DE results | `results_list_widget()` | `results_list_reactable()` |

Saving them:

```r
# A complete page: HTML plus a neighbouring lib/ directory. Keep them together.
page <- results_list_widget(res$DEPs_results, protein_quant = prep$protein_quant,
                            comparisons = "B-A")
htmltools::save_html(page, "results_B-A.html")

# A single self-contained file (widgets only: *_reactable() and the Highcharts)
tbl <- results_list_reactable(res$DEPs_results, comparisons = "B-A")
htmlwidgets::saveWidget(tbl, "results_B-A.html", selfcontained = TRUE)
```

The self-contained file travels alone but loses the NADIA toolbar. The Excel
export in `*_widget()` pages ships ExcelJS and PapaParse inside the package, so
it works offline as long as the `lib/` directory travels with the HTML.

## What to keep for reproducibility

- the original report and the sample metadata;
- the analysis script;
- a `.nadia` archive, or the exported tables;
- the processing parameters and software versions (`res$parameters`,
  `sessionInfo()`);
- any HTML dependency directory, if you shared interactive pages.

A `.nadia` file reduces the number of files to coordinate, but it complements
the original report and the script rather than replacing them.
