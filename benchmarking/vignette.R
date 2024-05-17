library(tidyverse)

taxonomy <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.tax"
fasta <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.fasta"

genera <- read_tsv(taxonomy,
                   col_names = c("accession", "taxonomy")) |>
  mutate(taxonomy = stringi::stri_replace_all_regex(taxonomy, ";$", ""))

fasta_data <- scan(fasta,
                   sep = "\n",
                   what = character(), quiet = TRUE)

sequence_names <- fasta_data[seq(1, length(fasta_data), by = 2)] |>
  str_replace_all(pattern = "^>([^\t]*)\t.*", "\\1")

sequences <- fasta_data[seq(2, length(fasta_data), by = 2)]

seq_table <- tibble(accession = sequence_names, sequence = sequences) |>
  inner_join(genera, by = "accession")

db <- build_kmer_database(seq_table$sequence,
                          seq_table$taxonomy,
                          kmer_size = 8)

unknown_sequence <- sequences[[1]]
bacteroides <- "TACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTTAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTATCAAACAGG"
oscillospiraceae <- "TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGCCGGGTTGACAAGTCAGATGTGAAATCCTGCGGCTTAACCGCAGAACTGCATTTGAAACTGTTGATCTTGAGTACTGGAGAGGCAGACGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGTCTGCTGGACAGCAACTGACGCTGAAGCACGAAAGTGCGGGGATCGAACAGG"
bacteroidales <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTTAAGTCAGTGGTCAAATTGAGGGGCTCAACCCCTTCCCGCCATTGAAACTGGCGATCTTGAGTGGAAGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATGCCGGCTTCCTACTGACGCTCATGCACGAAAGTGTGGGTAACGAACAGG"

num_bootstraps <- 1000
kmer_size <- 8

classify_sequence(unknown = unknown_sequence, database = db,
                  num_bootstraps = num_bootstraps, kmer_size = kmer_size)

profvis::profvis(
  classify_sequence(unknown = bacteroidales, database = db,
                    num_bootstraps = num_bootstraps, kmer_size = kmer_size)
)
#
filtered <- filter_taxonomy(consensus)


#
print_taxonomy(filtered)
