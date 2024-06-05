library(tidyverse)

fasta <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.fasta"
taxonomy <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.tax"


fasta_df <- read_fasta(fasta)

genera <- read_tsv(taxonomy,
                   col_names = c("accession", "taxonomy")) |>
  mutate(taxonomy = stringi::stri_replace_all_regex(taxonomy, ";$", ""))



seq_table <- as_tibble(fasta_df) |>
  inner_join(genera, by = c("id" = "accession"))

db <- build_kmer_database(seq_table$sequence,
                          seq_table$taxonomy,
                          kmer_size = 8)


unknown_sequence <- sequences[[1]]
bacteroides <- "TACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTTAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTATCAAACAGG"
oscillospiraceae <- "TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGCCGGGTTGACAAGTCAGATGTGAAATCCTGCGGCTTAACCGCAGAACTGCATTTGAAACTGTTGATCTTGAGTACTGGAGAGGCAGACGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGTCTGCTGGACAGCAACTGACGCTGAAGCACGAAAGTGCGGGGATCGAACAGG"
bacteroidales <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTTAAGTCAGTGGTCAAATTGAGGGGCTCAACCCCTTCCCGCCATTGAAACTGGCGATCTTGAGTGGAAGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATGCCGGCTTCCTACTGACGCTCATGCACGAAAGTGTGGGTAACGAACAGG"

num_bootstraps <- 100
kmer_size <- 8

consensus <- classify_sequence(unknown = bacteroidales, database = db,
                               num_bootstraps = num_bootstraps,
                               kmer_size = kmer_size)

filtered <- filter_taxonomy(consensus)
print_taxonomy(filtered)

