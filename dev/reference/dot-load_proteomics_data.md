# Create SummarizedExperiment from proteomics data

Create SummarizedExperiment from proteomics data

## Usage

``` r
.load_proteomics_data(
  data,
  metadata,
  protein_column = "ProteinGroups",
  gene_column = "GeneNames",
  condition_column = "Condition",
  label_column = "Column"
)
```

## Arguments

- data:

  Data frame with proteins and intensity values

- metadata:

  Data frame with sample information

- protein_column:

  Protein ID column name

- gene_column:

  Gene name column name

- condition_column:

  Condition column name in metadata

- label_column:

  Sample label column name in metadata

## Value

SummarizedExperiment with assays: raw, log2
