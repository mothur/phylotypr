# Classify a vector of 16S rRNA gene sequence fragments

The `classify_sequences()` function implements the Wang et al. naive
Bayesian classification algorithm for 16S rRNA gene sequence.

## Usage

``` r
classify_sequences(
  data = NULL,
  sequence,
  sequence_name = as.character(seq_along(sequence)),
  database,
  kmer_size = 8,
  num_bootstraps = 100,
  min_confidence = 80,
  seed = 1,
  num_processors = parallelly::availableCores()
)
```

## Arguments

- data:

  A data frame containing a column of `sequence_name` and `sequence`. If
  provided column names need to be provided to `sequence_name` and
  `sequence`.

- sequence:

  A vector of character objects representing DNA sequence that need to
  be classified. If `data` is provided then only the name of the column
  containing sequence needs to be provided, quotation marks are
  optional. This argument can be ignored when using a `{strollur}`
  object

- sequence_name:

  A vector of character objects containing the sequence name that need
  to be classified. Optional. Defaults to the position of the sequence
  in the vector of sequence. If `data` is provided then only the name of
  the column containing the sequence name needs to be provided,
  quotation marks are optional. This argument can be ignored when using
  a `{strollur}` object

- database:

  A kmer database generated using
  [`build_kmer_database`](https://mothur.org/phylotypr/reference/build_kmer_database.md)

- kmer_size:

  An integer value (default of 8) indicating the size of kmers to use
  for classifying sequence. Higher values use more RAM with potentially
  more specificity Lower values use less RAM with potentially less
  specificity. Benchmarking has found that the default of 8 provides the
  best specificity with the lowest possible memory requirement and
  fastest execution time.

- num_bootstraps:

  An integer value (default of 100). The value of `num_bootstraps` is
  the number of randomizations to perform where `1/kmer_size` of all
  kmers are sampled (without replacement) from `unknown_sequence`.
  Higher values will provide greater precision on the confidence score.

- min_confidence:

  An integer value (default of 80). The value of `min_confidence`
  represents the minimum number of bootstraps that must have a
  classification that agree with each other.

- seed:

  An integer value (default of 1). The value of `seed` sets the seed for
  the random number generator.

- num_processors:

  An integer value. The value of `num_processors` indicates the number
  of desired processors. Defaults to all available cores using
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)

## Value

A data frame object with three columns: the sequence name
(`sequence_name`), the taxonomic assignment for each level (`taxonomy`),
and the confidence associated with that classification level
(`confidence`). The value in `confidence` represents the percentage of
`num_bootstraps` that the classifier gave the same classification at
that level.

## References

Wang Q, Garrity GM, Tiedje JM, Cole JR. Naive Bayesian classifier for
rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl
Environ Microbiol. 2007 Aug;73(16):5261-7.
[doi:10.1128/AEM.00062-07](https://doi.org/10.1128/AEM.00062-07) PMID:
17586664; PMCID: PMC1950982.

## Examples

``` r
kmer_size <- 3
reference_names <- c("seq_1", "seq_2", "seq_3")
reference_sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
reference_genera <- c("A", "B", "B")
reference <- create_database(reference_names, reference_sequences, reference_genera)

db <- build_kmer_database(reference, kmer_size)

unknown_sequence <- c("ATGCGCTC", "ATGCGCTC", "ATGCGCTC")
unknown_id <- c("seq_1", "seq_2", "seq_3")

classify_sequences(
  sequence = unknown_sequence,
  database = db,
  kmer_size = kmer_size,
  seed = 1,
  num_processors = 1
)
#>   sequence_name level taxonomy confidence
#> 1             1     1        B        100
#> 2             2     1        B        100
#> 3             3     1        B        100

classify_sequences(
  sequence = unknown_sequence,
  sequence_name = unknown_id,
  database = db,
  kmer_size = kmer_size
)
#>   sequence_name level taxonomy confidence
#> 1         seq_1     1        B        100
#> 2         seq_2     1        B        100
#> 3         seq_3     1        B        100

classify_sequences(
  sequence = "ATGATCGAGTC",
  database = db,
  kmer_size = kmer_size,
  min_confidence = 95
)
#>   sequence_name level taxonomy confidence
#> 1             1     1        A         99

sequence_df <- data.frame(
  seq_names = unknown_id,
  seqs = unknown_sequence
)

classify_sequences(sequence_df,
  sequence_name = seq_names, sequence = seqs,
  database = db, kmer_size = 3
)
#>   sequence_name level taxonomy confidence
#> 1         seq_1     1        B        100
#> 2         seq_2     1        B        100
#> 3         seq_3     1        B        100
```
