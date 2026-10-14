# phylotypr

In the microbiome field, 16S rRNA gene sequencing has been a popular
method of characterizing the types of bacteria found in a community. The
most frequently used tool for classifying these sequences has been the
[naive Bayesian classifier](https://doi.org/10.1128/AEM.00062-07) that
used to be accessible through the Ribosomal Database Project and is
currently found in tools like mothur and QIIME. For the first time, the
[phylotypr](https://github.com/mothur/phylotypr) package provides an
R-based implementation of that tool. There are minimal dependencies and
the easy to use code executes rapidly.

## Getting started with phylotypr

To get started you’ll need to install
[phylotypr](https://github.com/mothur/phylotypr) using either of these
commands:

``` r

# devtools::install_github("mothur/phylotypr") #install development version
install.packages("phylotypr") # install stable version from CRAN
```

To get access to the packages function and a demonstration version of
the training data, you’ll need to run the
[`library()`](https://rdrr.io/r/base/library.html) function.

``` r

library(phylotypr)
```

## Building the database

Before we can classify a sequence, we need to generate the database.
This can be done using the
[`build_kmer_database()`](https://mothur.org/phylotypr/reference/build_kmer_database.md)
function and takes several seconds to execute.
[phylotypr](https://github.com/mothur/phylotypr) includes two databases:
`rdp_trainset9` and `rdp_trainset9_pds`. If you want others, you are
encouraged to install the `{phylotyprrefdata}` package from GitHub. You
can also create your own reference set using
[`create_database()`](https://mothur.org/phylotypr/reference/create_database.md).
For more details, see the “Alternative databases” section at the end of
this vignette. We’ll use `rdp_trainset9_pds` for this tutorial.

``` r

db <- build_kmer_database(rdp_trainset9_pds)
```

By default,
[`build_kmer_database()`](https://mothur.org/phylotypr/reference/build_kmer_database.md)
will use 8-nt long kmers to calculate its probabilities. This size can
be altered by the user if needed, but 8-nt performs very well and
shouldn’t need to be adjusted.

## Classifying a single unknown sequence

To classify a sequence we need the DNA sequence to be a character object
like `unknown` is below.

``` r

unknown <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTTAAGTCAGTGGTCAAATTGAGGGGCTCAACCCCTTCCCGCCATTGAAACTGGCGATCTTGAGTGGAAGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATGCCGGCTTCCTACTGACGCTCATGCACGAAAGTGTGGGTAACGAACAGG"
```

The `unknown` object can then be classified using
[`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md)
along with the database object from above (i.e., `db`). Because the
classification algorithm makes use of a random number generator, I
strongly encourage you to set the seed to get consistent results. This
needs to be set in
[`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md)
and not with [`set.seed()`](https://rdrr.io/r/base/Random.html). The
default seed is 1.

``` r

consensus <- classify_sequences(sequence = unknown, database = db, seed = 123)
```

The output of `classify_sequences` is a data frame with taxonomic levels
from kingdom (or domain) down to the genus and the confidence scores for
that sequence in each of those taxonomic levels

``` r

consensus
#>   sequence_name level      taxonomy confidence
#> 1             1     1      Bacteria        100
#> 2             1     2 Bacteroidetes         98
#> 3             1     3   Bacteroidia         82
#> 4             1     4 Bacteroidales         82
```

In this example, you’ll notice that the family and genus assignment have
confidence scores below 80%. This is the value commonly used in mothur
and which was used on the RDP website. To filter the classification to
that level of confidence, you can use the
[`filter_taxonomy()`](https://mothur.org/phylotypr/reference/filter_taxonomy.md)
function:

``` r

filtered <- filter_taxonomy(consensus)
```

You should see that the “Porphyromonadaceae” and “Barnesiella” names are
removed along with its confidence score.

``` r

filtered
#>   sequence_name level      taxonomy confidence
#> 1             1     1      Bacteria        100
#> 2             1     2 Bacteroidetes         98
#> 3             1     3   Bacteroidia         82
#> 4             1     4 Bacteroidales         82
```

To convert that list format to a string, you can use the
[`print_taxonomy()`](https://mothur.org/phylotypr/reference/print_taxonomy.md)
function:

``` r

print_taxonomy(filtered)
#>   sequence_name
#> 1             1
#>                                                            taxonomy
#> 1 Bacteria(100);Bacteroidetes(98);Bacteroidia(82);Bacteroidales(82)
```

## Classifying multiple unknown sequences

Most likely you will have multiple sequences that you would like to
classify. The [phylotypr](https://github.com/mothur/phylotypr) package
installs with an example FASTA-formatted file, “miseq_sop.fasta”. The
code below shows you can read in a FASTA-formatted file and classify the
sequences. Note that for classifying your own data you would replace
`phylotypr_example(...)` with the path to your own data. By default
[`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md)
will use as many processors as your computer will allow. We will show
how to limit the number of processors later in this vignette:

``` r

read_fasta(phylotypr_example("miseq_sop.fasta.gz")) |>
  classify_sequences(database = db, seed = 123) |>
  print_taxonomy() |>
  head()
```

    #>                                 sequence_name
    #> 1 M00967_43_000000000A3JHG_1_1101_10331_23332
    #> 2 M00967_43_000000000A3JHG_1_1101_10382_22128
    #> 3 M00967_43_000000000A3JHG_1_1101_10491_15766
    #> 4  M00967_43_000000000A3JHG_1_1101_10877_9025
    #> 5 M00967_43_000000000A3JHG_1_1101_11035_15765
    #> 6  M00967_43_000000000A3JHG_1_1101_11296_4480
    #>                                                                                                                        taxonomy
    #> 1            Bacteria(100);Firmicutes(98);Clostridia(98);Clostridiales(98);Lachnospiraceae(88);Lachnospiraceae_unclassified(88)
    #> 2         Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(99);Ruminococcaceae_unclassified(99)
    #> 3 Bacteria(100);Bacteroidetes(100);Bacteroidia(99);Bacteroidales(99);Porphyromonadaceae(96);Porphyromonadaceae_unclassified(96)
    #> 4         Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100)
    #> 5         Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(97);Lachnospiraceae_unclassified(97)
    #> 6            Bacteria(100);Firmicutes(95);Clostridia(94);Clostridiales(93);Ruminococcaceae(88);Ruminococcaceae_unclassified(88)

Alternatively, you may have your data stored as a
[strollur](https://github.com/mothur/strollur) object you can do the
following

``` r

strollur::new_dataset(dataset_name = "sop") |>
  strollur::add(
    table = phylotypr::read_fasta(phylotypr_example("miseq_sop.fasta.gz"))
  ) |>
  classify_sequences(db, seed = 123) |>
  print_taxonomy() |>
  head()
```

    #> Added 2490 sequences.
    #> Assigned 2490 sequence taxonomies.
    #>                                 sequence_name
    #> 1 M00967_43_000000000A3JHG_1_1101_10331_23332
    #> 2 M00967_43_000000000A3JHG_1_1101_10382_22128
    #> 3 M00967_43_000000000A3JHG_1_1101_10491_15766
    #> 4  M00967_43_000000000A3JHG_1_1101_10877_9025
    #> 5 M00967_43_000000000A3JHG_1_1101_11035_15765
    #> 6  M00967_43_000000000A3JHG_1_1101_11296_4480
    #>                                                                                                                        taxonomy
    #> 1            Bacteria(100);Firmicutes(98);Clostridia(98);Clostridiales(98);Lachnospiraceae(88);Lachnospiraceae_unclassified(88)
    #> 2         Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(99);Ruminococcaceae_unclassified(99)
    #> 3 Bacteria(100);Bacteroidetes(100);Bacteroidia(99);Bacteroidales(99);Porphyromonadaceae(96);Porphyromonadaceae_unclassified(96)
    #> 4         Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100)
    #> 5         Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(97);Lachnospiraceae_unclassified(97)
    #> 6            Bacteria(100);Firmicutes(95);Clostridia(94);Clostridiales(93);Ruminococcaceae(88);Ruminococcaceae_unclassified(88)

## Classifying multiple unknown sequences in parallel

The [phylotypr](https://github.com/mothur/phylotypr) package allows you
to run
[`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md)
in parallel. So, if your computer has 8 processors,
[`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md)
will split your data into 8 chunks and then classify each chunk
separately before pulling it all together. To see how many processors
you have available, you can run the following:

``` r

parallelly::availableCores()
#> system 
#>      4
```

By default,
[`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md)
will use all of those processors. If you need to limit the number of
processors, you can use the `num_processors` argument

``` r

read_fasta(phylotypr_example("miseq_sop.fasta.gz")) |>
  classify_sequences(database = db, seed = 123, num_processors = 2) |>
  print_taxonomy() |>
  head()
#>                                 sequence_name
#> 1 M00967_43_000000000A3JHG_1_1101_10331_23332
#> 2 M00967_43_000000000A3JHG_1_1101_10382_22128
#> 3 M00967_43_000000000A3JHG_1_1101_10491_15766
#> 4  M00967_43_000000000A3JHG_1_1101_10877_9025
#> 5 M00967_43_000000000A3JHG_1_1101_11035_15765
#> 6  M00967_43_000000000A3JHG_1_1101_11296_4480
#>                                                                                                                        taxonomy
#> 1            Bacteria(100);Firmicutes(98);Clostridia(98);Clostridiales(98);Lachnospiraceae(88);Lachnospiraceae_unclassified(88)
#> 2         Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(99);Ruminococcaceae_unclassified(99)
#> 3 Bacteria(100);Bacteroidetes(100);Bacteroidia(99);Bacteroidales(99);Porphyromonadaceae(96);Porphyromonadaceae_unclassified(96)
#> 4         Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100)
#> 5         Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(97);Lachnospiraceae_unclassified(97)
#> 6            Bacteria(100);Firmicutes(95);Clostridia(94);Clostridiales(93);Ruminococcaceae(88);Ruminococcaceae_unclassified(88)
```

## Run-to-run variation

Consider the results of classifying `unknown` three times with three
different seeds:

``` r

classify_sequences(sequence = unknown, database = db, seed = 1) |>
  print_taxonomy()
#>   sequence_name
#> 1             1
#>                                                                                   taxonomy
#> 1 Bacteria(100);Bacteroidetes(96);Bacteroidia(85);Bacteroidales(85);Porphyromonadaceae(84)
classify_sequences(sequence = unknown, database = db, seed = 2) |>
  print_taxonomy()
#>   sequence_name
#> 1             1
#>                                                                                   taxonomy
#> 1 Bacteria(100);Bacteroidetes(95);Bacteroidia(83);Bacteroidales(83);Porphyromonadaceae(82)
classify_sequences(sequence = unknown, database = db, seed = 3) |>
  print_taxonomy()
#>   sequence_name
#> 1             1
#>                                                            taxonomy
#> 1 Bacteria(100);Bacteroidetes(92);Bacteroidia(80);Bacteroidales(80)
```

You’ll notice that the classification of the three replicates varies
slightly. Again, this is due to the use of the random number generator
and getting confidence scores that bounce around 80%, which results in
different levels of filtering. To get a more stable classification, you
should be setting your random number generator seed with a consistent
value of `seed`.

To get a more precise confidence score, you could set
`num_bootstraps = 1000` in
[`classify_sequences()`](https://mothur.org/phylotypr/reference/classify_sequences.md).
Be forewarned that the function will take 10-times longer to run with
10-times the number of bootstraps!

``` r

classify_sequences(sequence = unknown, database = db, num_bootstraps = 1000,
                   seed = 123) |>
  print_taxonomy()
#>   sequence_name
#> 1             1
#>                                                                   taxonomy
#> 1 Bacteria(99.9);Bacteroidetes(96.5);Bacteroidia(83.9);Bacteroidales(83.9)
```

## Alternative databases

The [phylotypr](https://github.com/mothur/phylotypr) package ships with
the RDP’s v.9 of their training data. This is relatively small and old
(2010) relative to their latest versions. You are encouraged to install
newer versions of the RDP, greengenes, and SILVA databases from the
`{phylotyprrefdata}` package on GitHub. Note that installing the package
will take about 20 minutes to install. If it sits at “moving datasets to
lazyload DB” for a long time, this is normal :)

``` r

devtools::install_github("mothur/phylotyprrefdata")
library(phylotyprrefdata)
```

The following will list the references that are available in
`{phylotyprrefdata}`:

``` r

data(package = "phylotyprrefdata")
```
