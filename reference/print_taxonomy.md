# Print taxonomy for an unknown sequence

The `print_taxonomy()` will output the consensus taxonomy for an unknown
sequence with confidence scores for each taxonomic level and each
taxonomic level separated by semi-colons

## Usage

``` r
print_taxonomy(
  consensus,
  n_levels = 0,
  fill_levels = TRUE,
  show_confidence = TRUE
)
```

## Arguments

- consensus:

  A data frame that contains at least three columns: the `sequence_name`
  column that contains an identifier for the sequence, the `taxonomy`
  column contains the classification at each taxonomic level, and the
  `confidence` column contains the percentage of bootstraps that had the
  specified classification. The `level` column is optional and contains
  the numeric level of the level (e.g., Kingdom = 1, Genus = 6)

- n_levels:

  An integer indicating the number of taxonomic levels to expect. If the
  number of observed levels is less than this value, and `fill_levels`
  is `TRUE` then missing levels will have "\_unclassified" appended to
  the end of the last named classification until the desired level is
  reached. If `n_levels` is less than the observed number of levels,
  then the string of the output will be truncated to `n_levels`. If 0 is
  given then we will report the longest possible taxonomy. Default is 0.

- fill_levels:

  A logical indicating whether or not to extend the classifications of
  sequences that don't meet `n_levels`. Default is `TRUE`.

- show_confidence:

  A logical indicating whether or not to include the confidence scores.
  Default is `TRUE`

## Value

A vector of character strings indicating the classification at each
taxonomic level with the corresponding confidence in parentheses (if
`show_confidence` is `FALSE`. Each taxonomic level is separated by a
semi-colon

## Examples

``` r
oscillospiraceae <- data.frame(
  sequence_name = c("seq_1", "seq_1", "seq_1", "seq_1", "seq_1"),
  level = 1:5,
  taxonomy = c(
    "Bacteria", "Bacillota", "Clostridia",
    "Eubacteriales", "Oscillospiraceae"
  ),
  confidence = c(100, 100, 99, 99, 98)
)

print_taxonomy(oscillospiraceae, n_levels = 6)
#>   sequence_name
#> 1         seq_1
#>                                                                                                               taxonomy
#> 1 Bacteria(100);Bacillota(100);Clostridia(99);Eubacteriales(99);Oscillospiraceae(98);Oscillospiraceae_unclassified(98)
print_taxonomy(oscillospiraceae, n_levels = 4)
#>   sequence_name                                                      taxonomy
#> 1         seq_1 Bacteria(100);Bacillota(100);Clostridia(99);Eubacteriales(99)
print_taxonomy(oscillospiraceae, n_levels = 6, show_confidence = FALSE)
#>   sequence_name
#> 1         seq_1
#>                                                                                     taxonomy
#> 1 Bacteria;Bacillota;Clostridia;Eubacteriales;Oscillospiraceae;Oscillospiraceae_unclassified
print_taxonomy(oscillospiraceae, n_levels = 7, fill_levels = FALSE)
#>   sequence_name
#> 1         seq_1
#>                                                                             taxonomy
#> 1 Bacteria(100);Bacillota(100);Clostridia(99);Eubacteriales(99);Oscillospiraceae(98)
```
