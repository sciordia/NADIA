# Import Protein Intensities from a SummarizedExperiment

Converts a `SummarizedExperiment` of protein intensities into the
`proteomics_data` object that every NADIA reader returns. It is the way
in for data imported with other Bioconductor tools, for instance one of
the assays of a `QFeatures` object.
[`process_proteomics()`](https://sciordia.github.io/NADIA/dev/reference/process_proteomics.md)
also accepts the `SummarizedExperiment` directly and calls this
function; call it yourself to keep the object, for example to store it
with
[`write_nadia()`](https://sciordia.github.io/NADIA/dev/reference/write_nadia.md).

## Usage

``` r
preprocess_se(
  se,
  condition_order = NULL,
  assay_name = NULL,
  log2_input = FALSE,
  condition_column = "Condition",
  replicate_column = "Replicate",
  verbose = TRUE
)
```

## Arguments

- se:

  A `SummarizedExperiment` with protein intensities.

- condition_order:

  Optional character vector with the conditions to keep, in the order
  they should appear. `NULL` (default) keeps them all, in the order of
  the factor levels or, for text, of first appearance.

- assay_name:

  Name or index of the assay with the intensities. `NULL` (default)
  takes the first one.

- log2_input:

  `TRUE` if the intensities are on the log2 scale. Default `FALSE`,
  linear intensities.

- condition_column:

  Column of `colData(se)` with the condition of each sample. Default
  `"Condition"`.

- replicate_column:

  Column of `colData(se)` with the replicate number. Default
  `"Replicate"`; if the column is absent, the replicates are numbered.

- verbose:

  Print progress messages. Default `TRUE`.

## Value

An object of class `c("se_data", "proteomics_data", "list")` with the
elements `metadata`, `protein_id` and `protein_quant`, like the other
readers.

## Details

What is read from the object:

- **Intensities**: the assay `assay_name`, proteins in rows and samples
  in columns. With `log2_input = TRUE` the values are taken back to the
  linear scale, because NADIA applies the log2 transformation when it
  normalises. Zeros are treated as missing values, as for the other
  readers.

- **Proteins**: the row names, which must be present and unique. The
  gene names come from `rowData()$Genes` and the descriptions from
  `rowData()$Description`, when those columns exist.

- **Design**: the condition from `colData()[[condition_column]]` and the
  replicate from `colData()[[replicate_column]]`; without a replicate
  column the samples are numbered within each condition, in column
  order.

The samples are renamed `<condition>_<replicate>`, the convention the
rest of NADIA relies on, and the original column names are kept in
`metadata$R.FileName`. Identification metrics (precursor or peptide
counts, scores) are not read from `rowData()`, so those columns of the
tables are left empty.

## See also

[`process_proteomics()`](https://sciordia.github.io/NADIA/dev/reference/process_proteomics.md),
which accepts the `SummarizedExperiment` directly.

## Examples

``` r
data(nadia_dia)
# A SummarizedExperiment built from the example data
quant <- as.matrix(nadia_dia$protein_quant[, grep("^PG\\.Quantity_",
                                                  names(nadia_dia$protein_quant))])
rownames(quant) <- nadia_dia$protein_quant$PG.ProteinGroups
se <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(intensity = quant),
  colData = data.frame(Condition = nadia_dia$metadata$R.Condition,
                       Replicate = nadia_dia$metadata$R.Replicate),
  rowData = data.frame(Genes = nadia_dia$protein_quant$PG.Genes))

pre <- preprocess_se(se, verbose = FALSE)
pre
#> Protein intensities imported from a SummarizedExperiment
#> --------------------------------------------------------
#> Samples (metadata): 12 
#> Proteins: 2000 
#> 
#> Conditions: A, B, D 

# The same values on the log2 scale
SummarizedExperiment::assay(se, "intensity") <- log2(quant)
pre_log2 <- preprocess_se(se, log2_input = TRUE, verbose = FALSE)
```
