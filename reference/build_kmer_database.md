# Build kmer database

Build kmer database for classifying 16S rRNA and other gene sequences to
a genus when a kmer size is provided.

## Usage

``` r
build_kmer_database(classification_database, kmer_size = 8)
```

## Arguments

- classification_database:

  A `classification_database` object, which is a list containing the
  fields `seq_taxonomy_data` and `reference_data`. These are both
  data.frame objects. The package includes rdp_trainset9 and
  rdp_trainset9_pds. Other databases can be obtained from
  phylotyprrefdata or by using
  [`create_database()`](https://mothur.org/phylotypr/reference/create_database.md)

- kmer_size:

  An integer indicating the length of the nucleotide word to base our
  classification on (default = 8)

## Value

A list object containing the genus level conditional probability
(`conditional_prob`) of seeing each kmer in a given genus as well as the
genus names (`genera`)

## References

Wang Q, Garrity GM, Tiedje JM, Cole JR. Naive Bayesian classifier for
rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl
Environ Microbiol. 2007 Aug;73(16):5261-7.
[doi:10.1128/AEM.00062-07](https://doi.org/10.1128/AEM.00062-07) PMID:
17586664; PMCID: PMC1950982.

## Examples

``` r
kmer_size <- 8

if (FALSE) { # \dontrun{
# Function takes >5 s to run on this database
build_kmer_database(rdp_trainset9_pds, kmer_size)
} # }
```
