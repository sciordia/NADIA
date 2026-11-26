# Filter proteins by group presence

Keeps proteins that have enough non-NA values in at least a minimum
number of experimental groups.

## Usage

``` r
.filter_proteins_by_group(
  data,
  metadata,
  min_reps = NULL,
  min_groups = 1,
  grouping_column = "Condition"
)
```

## Arguments

- data:

  Intensity matrix (proteins x samples)

- metadata:

  Data frame with sample information

- min_reps:

  Minimum replicates with non-NA values per group. If NULL, uses half of
  the smallest group size.

- min_groups:

  Minimum groups meeting min_reps (default: 1)

- grouping_column:

  Name of grouping column in metadata (default: "Condition")

## Value

List with:

- data: Filtered matrix

- keep: Logical vector of kept rows

- summary: Filtering summary
