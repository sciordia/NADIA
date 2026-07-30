# NADIA

**Missing Value-Aware DIA Proteomics Analysis** — an R package for differential
protein expression analysis, with particular attention to the missing values
(**NA**) that characterise data-independent acquisition (**DIA**).

## Why the missing values

A DIA experiment routinely leaves 5–20 % of its quantification matrix empty, and
those blanks are not an accident of the instrument. A protein below the detection
limit in one condition and above it in another produces exactly the pattern a
differential expression test is supposed to find. Deleting those proteins throws
away the signal; filling them in carelessly manufactures it.

The structure is visible in the example dataset shipped with the package. Cross
the percentage of missing values against the differential expression call:

| missing | Up | Down | No Change |
|---|---|---|---|
| 0 % | 795 | 209 | 3761 |
| 1–25 % | 90 | 57 | 239 |
| **26–50 %** | **70** | **619** | 54 |
| > 50 % | 30 | 52 | 15 |

Proteins with no gaps are mostly unchanged. Proteins missing between a quarter
and a half of their values are called *down* nine times out of ten. That
asymmetry is not a bug in the test — it is the signature of values **missing not
at random**, and it is what the package is built to handle rather than ignore.

## What it covers

From the raw report to the interactive figure:

1. **Preprocessing** — Spectronaut/DIA-NN, TMT (Proteome Discoverer) and LFQ
   (Proteome Discoverer); all three produce the same `proteomics_data` S3 object,
   so the rest of the pipeline is indifferent to the source.
2. **Processing** — normalisation (13 methods), optional batch correction (BERT,
   applying ComBat, limma or a reference batch, plus PVCA diagnostics), imputation
   (19 methods, including the MAR/MNAR hybrids `combo` and `softHybrid` and the
   probabilistic model `limpa`) and differential expression (`limma` or `limpa`).
3. **Metrics and benchmarking** — normalisation assessment (PCV/PMAD/PEV,
   intragroup correlation, group separation) and imputation assessment (NAguideR
   framework: NRMSE, SOR, PSS, ACC_OI), plus benchmarking on *spike-in* datasets
   and OpDEA ranking of normalisation × imputation combinations.
4. **Visualisation** — interactive with Highcharts (boxplots, volcano, PCA,
   cluster profiles) and static with ggplot2/ComplexHeatmap, along with
   interactive tables built on reactable.

## Method choice is not cosmetic

Four pipelines, same data, same test, same threshold:

| pipeline | Up | Down | No Change |
|---|---|---|---|
| `cycloess` + `combo` | 985 | 937 | 4069 |
| `cycloess` + `knn` | 872 | 179 | 4940 |
| `quantile` + `combo` | 902 | 2825 | 2264 |
| `log2` + `min` | 761 | 778 | 4452 |

The count of down-regulated proteins spans a factor of sixteen. That is the
argument for measuring the choice rather than inheriting it — and for stating in
any publication which normalisation and imputation produced the table, since
without that the numbers are not reproducible even from the same raw data.

`vignette("choosing-methods")` shows how to score the alternatives when you have
no ground truth; `vignette("benchmarking")` when a spike-in gives you one.

## Installation

NADIA is an R package and requires R >= 4.4. It is not on Bioconductor yet, so
install it from this repository:

```r
# install.packages("remotes")
remotes::install_github("sciordia/NADIA", build_vignettes = TRUE)
library(NADIA)
```

The required dependencies (the `Imports:` field) are installed automatically. The
**optional** ones (`Suggests:`) are only needed if you use the method that calls
them — `mice` only with `imp_method = "mice"`, `pROC` for the AUC/pAUC
benchmarking metrics, `Mfuzz` for the Pattern Profiler. When one is missing, the
function says which one. To install them all at once:

```r
source("install_dependencies.R")         # from a clone of the repository
install_nadia_deps(dry_run = TRUE)       # only report what is missing
install_nadia_deps(optional = FALSE)     # only the essentials
```

## A first example

```r
library(NADIA)

# Preprocessed example dataset: 3 conditions x 4 replicates, 2,000 proteins,
# 8 % missing values (a random sample of the full experiment, so the missingness
# and the differential-expression balance are representative)
data(nadia_dia)

res <- process_proteomics(nadia_dia,
                          norm_method = "cycloess",
                          imp_method  = "combo",
                          de_method   = "limma")
head(res$DEPs_results)

# Or starting from the raw report
prep <- preprocess_spectronaut(
  system.file("extdata", "nadia_dia_report.tsv.gz", package = "NADIA"),
  condition_order = c("A", "B", "D"))
```

`process_proteomics()` writes nothing to disk unless it is given an `export_dir`.

Alongside `logFC` and `adj.P.Val`, the results carry `MissingGlobal`,
`MissingPCT1` and `MissingPCT2` — the percentage of missing values in each group
**before** imputation. They are what lets you tell a real on/off signal from a
fold change built entirely on imputed values, which no p-value can do on its own.

## Documentation

Nine vignettes cover the pipeline end to end. Start with `vignette("NADIA")`.

| Vignette | Question it answers |
|---|---|
| **`NADIA`** | How do I analyse my data? |
| `input-formats` | My data are TMT or label-free, not DIA |
| `missing-values` | Which imputation method, and why |
| `choosing-methods` | How do I pick, with no ground truth |
| `benchmarking` | I have a spike-in, so I *do* have ground truth |
| `batch-correction` | My samples were run in batches |
| `pattern-profiler` | Which proteins behave alike across conditions |
| `visualization` | I want to adjust the figures |
| `results-and-export` | How do I get results out of R |

All of them execute real code on the shipped dataset; none is a static document.

## Package structure

```
NADIA/
├── R/                  22 files, 102 exported functions
├── man/                generated with roxygen2 -- do not edit by hand
├── NAMESPACE           generated with roxygen2 -- do not edit by hand
├── vignettes/          the nine articles listed above
├── tests/testthat/     453 assertions in 7 files
├── data/               nadia_dia, the preprocessed example dataset
├── inst/
│   ├── extdata/        trimmed DIA, TMT and LFQ reports
│   ├── scripts/        make_extdata.R -- how those reports were produced
│   ├── js/             ExcelJS + PapaParse (MIT), for the offline Excel export
│   └── css/            styles for the reactable tables
├── data-raw/           full datasets .......... not in the package
├── results/            real analysis outputs .. not in the package
├── CLAUDE.md           architecture, module by module
└── CODE_REVIEW_*.md    code reviews and the measured impact of their fixes
```

`data-raw/`, `results/`, the numbered scripts and the Markdown reports are
excluded from the build via `.Rbuildignore`.

### The modules

The code follows the shape of the pipeline. Each stage is one file with one main
exported function, so a stage can be run on its own or swapped out.

**Preprocessing** — raw report to a common S3 object

| File | Entry point |
|---|---|
| `Preprocessing.R` | `preprocess_spectronaut()` — Spectronaut/DIA-NN, long format |
| `Preprocessing_TMT.R` | `preprocess_tmt()` — design read from the `Abundance:` column suffixes |
| `Preprocessing_LFQ.R` | `preprocess_lfq()` — design read from a separate annotation file |

All three return `c("<type>_data", "proteomics_data", "list")` with the same three
elements (`metadata`, `protein_id`, `protein_quant`).

**Processing** — the core

| File | Entry point |
|---|---|
| `Processing.R` | `process_proteomics()` — the coordinator |
| `Normalization.R` | `normalize_proteomics()` — 13 methods |
| `Batch_Correction.R` | `batch_correct_proteomics()`, `pvca_analysis()` |
| `Imputation.R` | `impute_proteomics()` — 19 methods |
| `DEAnalysis.R` | `de_analysis_proteomics()` — limma or limpa |

`process_proteomics()` runs them in order — normalisation → optional batch
correction → imputation → differential expression — and returns the
`SummarizedExperiment` with one assay per stage, plus the results table and the
inputs the plotting modules expect.

**Metrics and benchmarking** — scoring the choices

| File | Entry point |
|---|---|
| `Normalization_Metrics.R` | `normalization_metrics()` — no ground truth needed |
| `Imputation_Metrics.R` | `imputation_metrics()` — ground truth by simulation |
| `Benchmarking_Single.R` | `benchmarking_proteomics()` — one pipeline against a spike-in |
| `Benchmarking_Multiple.R` | `benchmarking_multiple()` — OpDEA ranking across pipelines |

**Visualisation** — independent of the pipeline, they take data frames

| File | Entry point |
|---|---|
| `Volcano_Plot_Highcharts_Final.R` | `volcano_highchart_list()` |
| `Boxplot_Highcharts_Final.R` | `boxplot_highchart_list()` |
| `PCA_Highcharts_Final.R` | `pca_highchart_list()` |
| `Heatmap_tidyHeatmap.R` | `proteomics_heatmap()` — static, ComplexHeatmap |
| `Pattern_Profiler_Analysis.R` | `pattern_profiler_analysis()` — Mfuzz soft clustering |
| `Pattern_Profiler_Highcharts.R` | `cluster_profile_highchart_list()` |
| `Results_List_reactable.R` | `results_list_widget()` and the other three tables |

**Infrastructure**

| File | Contents |
|---|---|
| `NADIA-package.R` | every `@importFrom`, and the NSE column names |
| `utils.R` | shared internals: `%||%`, RNG, colour and filtering helpers |
| `data.R` | documentation for `nadia_dia` |

## Data

The example data come from quantitative proteomics experiments acquired at the
Proteomics Facility of the Centro Nacional de Biotecnologia (CNB-CSIC). The
reports in `inst/extdata/` are trimmed to a random sample of 2,000 protein groups
— random rather than best-covered on purpose, so that the missingness the package
exists to handle survives into the examples.

## License

GPL-3 (or any later version) (c) 2025 Sergio Ciordia
