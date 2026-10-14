# Build custom database object for use with phylotypr

The `create_database()` function builds a `classification_database`
object that can be used to generate the kmer database for downstream use
in
[`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md).
This function is only necessary if you are building your own database or
modifying one of the existing databases provided in `{phylotypr}` or
`{phylotyprrefdata}`.

## Usage

``` r
create_database(sequence_name, sequence, taxonomy, reference_data = NULL)
```

## Arguments

- sequence_name:

  A character vector of identifiers for each sequence in the database.
  Should have the same number of values as found in the `sequence` and
  `taxonomy` vectors.

- sequence:

  A character vector of the sequences that are included in the database.
  Should have the same number of values as found in the `sequence_name`
  and `taxonomy` vectors.

- taxonomy:

  A character vector of taxonomy data for each sequence in the database.
  Should have the same number of values as found in the `sequence_name`
  and `sequence` vectors. The format of each value should be with each
  taxonomic level separated by a semicolon. Ideally each value in
  `taxonomy` will have the same number of levels.

- reference_data:

  single row data frame object that has columns vendor, name, version,
  usage, notes, method_url, documentation_url, and citation. This is
  supplied with databases obtained from phylotyprrefdata. This
  information is only used if user supplied sequences are supplied as a
  strollur object.

## Value

A `classification_database` object, which is a list containing the
fields `seq_taxonomy_data` and `reference_data`. These are both
data.frame objects.

## Examples

``` r
sequence_name <- c("seq_a", "seq_b", "seq_c")
sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

create_database(
  sequence_name = sequence_name,
  sequence = sequence,
  taxonomy = genera
)
#> $seq_taxonomy_data
#>   sequence_name    sequence taxonomy.1 taxonomy.2 taxonomy.3 taxonomy.4
#> 1         seq_a AAAATCGAGTA          A          b          C          d
#> 2         seq_b ATGCGCTCCCC          A          b          C          e
#> 3         seq_c ATGCGCTCCCC          A          b          C          e
#> 
#> $reference_data
#> NULL
#> 
#> attr(,"class")
#> [1] "classification_database"
```
