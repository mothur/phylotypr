# Filter taxonomy

The `filter_taxonomy()` function will filter a consensus taxonomy to
remove any taxonomic levels where the confidence score is below a
`min_confidence` level

## Usage

``` r
filter_taxonomy(consensus, min_confidence = 80)
```

## Arguments

- consensus:

  A data frame that contains at least three columns: the `sequence_name`
  column that contains an identifier for the sequence, the `taxonomy`
  column contains the classification at each taxonomic level, and the
  `confidence` column contains the percentage of bootstraps that had the
  specified classification. The `level` column is optional and contains
  the numeric level of the level (e.g., Kingdom = 1, Genus = 6)

- min_confidence:

  A double value between 0 and 100 (default = 80). The minimum
  percentage of bootstrap replicates that had the same classification.
  Any confidence score below this value will have the corresponding
  taxonomy removed

## Value

A data frame object that is filtered to remove low confidence
taxonomies. One column, `taxonomy`, contains the taxonomy at each
taxonomic level and the other column, `confidence` contains the
confidence score for that taxonomy. There will be no taxonomies or
confidence scores below `min_confidence`

## Examples

``` r
oscillospiraceae <- data.frame(
  sequence_name = c("seq_1", "seq_1", "seq_1", "seq_1", "seq_1", "seq_1"),
  level = 1:6,
  taxonomy = c(
    "Bacteria", "Bacillota", "Clostridia",
    "Eubacteriales", "Oscillospiraceae",
    "Flintibacter"
  ),
  confidence = c(100, 100, 99, 99, 98, 58)
)

filter_taxonomy(oscillospiraceae, min_confidence = 80)
#>   sequence_name level         taxonomy confidence
#> 1         seq_1     1         Bacteria        100
#> 2         seq_1     2        Bacillota        100
#> 3         seq_1     3       Clostridia         99
#> 4         seq_1     4    Eubacteriales         99
#> 5         seq_1     5 Oscillospiraceae         98
```
