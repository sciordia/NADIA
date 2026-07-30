# NADIA

**Missing Value-Aware DIA Proteomics Analysis** — an R package for differential
protein expression analysis, with particular attention to the missing values
(**NA**) that characterise data-independent acquisition (**DIA**).

## What it covers

From the raw report to the interactive figure:

1. **Preprocessing** — Spectronaut/DIA-NN, TMT (Proteome Discoverer) and LFQ
   (Proteome Discoverer); all three produce the same `proteomics_data` S3 object.
2. **Processing** — normalization (13 methods), optional batch correction (BERT,
   applying ComBat, limma or a reference batch, plus PVCA diagnostics), imputation
   (19 methods, including the MAR/MNAR hybrids `combo` and `softHybrid` and the
   probabilistic model `limpa`) and differential expression (`limma` or `limpa`).
3. **Metrics and benchmarking** — normalization assessment (PCV/PMAD/PEV,
   intragroup correlation, group separation) and imputation assessment (NAguideR
   framework: NRMSE, SOR, PSS, ACC_OI), plus benchmarking on *spike-in* datasets
   and OpDEA ranking of normalization x imputation combinations.
4. **Visualization** — interactive with Highcharts (boxplots, volcano, PCA,
   cluster profiles) and static with ggplot2/ComplexHeatmap, along with
   interactive tables built on reactable.

## Installation

NADIA is an R package and requires R >= 4.4. It is not on Bioconductor yet, so
install it from this repository:

```r
# install.packages("remotes")
remotes::install_github("sciordia/NADIA")
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

`renv` is not used.

## A first example

```r
library(NADIA)

# Preprocessed example dataset: 3 conditions x 4 replicates, 2,000 proteins
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

## Repository layout

- `R/` — package code: 22 files, 102 exported functions.
- `man/`, `NAMESPACE` — generated with roxygen2; do not edit by hand.
- `inst/extdata/` — trimmed Spectronaut, TMT and LFQ reports for the examples;
  `inst/scripts/make_extdata.R` documents how they were obtained.
- `data/` — the example dataset `nadia_dia`.
- `data-raw/`, `results/` — full datasets and real analysis outputs. Not part of
  the package (`.Rbuildignore`).
- `example_workflow.R` and the numbered scripts — end-to-end walkthroughs over the
  full datasets.
- `CLAUDE.md` — detailed description of the architecture and of every module.
- `CODE_REVIEW_*.md` — code reviews and impact analyses of their fixes.

## Data

The example data come from quantitative proteomics experiments acquired at the
Proteomics Facility of the Centro Nacional de Biotecnologia (CNB-CSIC).

## License

MIT (c) 2025 Sergio Ciordia
