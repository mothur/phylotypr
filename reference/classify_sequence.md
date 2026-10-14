# Classify 16S rRNA gene sequence fragment

The `classify_sequence()` function implements the Wang et al. naive
Bayesian classification algorithm for 16S rRNA gene sequences.

## Usage

``` r
classify_sequence(sequence, database, kmer_size = 8, num_bootstraps = 100)
```

## Arguments

- sequence:

  A character object representing a DNA sequence that needs to be
  classified.

- database:

  A kmer database generated using
  [`build_kmer_database`](https://mothur.org/phylotypr/reference/build_kmer_database.md)

- kmer_size:

  An integer value (default of 8) indicating the size of kmers to use
  for classifying sequences. Higher values use more RAM with potentially
  more specificity Lower values use less RAM with potentially less
  specificity. Benchmarking has found that the default of 8 provides the
  best specificity with the lowest possible memory requirement and
  fastest execution time.

- num_bootstraps:

  An integer value (default of 100). The value of `num_bootstraps` is
  the number of randomizations to perform where `1/kmer_size` of all
  kmers are sampled (without replacement) from `sequence`. Higher values
  will provide greater precision on the confidence score.

## Value

A list object of two vectors. One vector (`taxonomy`) is the taxonomic
assignment for each level. The second vector (`confidence`) is the
percentage of `num_bootstraps` that the classifier gave the same
classification at that level.

## Examples

``` r
kmer_size <- 3
reference_names <- c("seq_1", "seq_2", "seq_3")
reference_sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
reference_genera <- c("A", "B", "B")
reference <- create_database(reference_names, reference_sequences, reference_genera)

db <- build_kmer_database(reference, kmer_size)
unknown_sequence <- "ATGCGCTC"

classify_sequence(
  sequence = unknown_sequence,
  database = db,
  kmer_size = kmer_size
)
#> Warning: 'classify_sequence' is deprecated.
#> Use 'classify_sequences' instead.
#> See help("Deprecated")
#>   sequence_name level taxonomy confidence
#> 1             1     1        B        100
```
