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
num_bootstraps <- 100
kmer_size <- 8

#classify_sequence(unknown = unknown_sequence, database = db)
kmers <- detect_kmers(sequence = unknown_sequence, kmer_size)

bs_class <- numeric(length = num_bootstraps)

for(i in 1:num_bootstraps){
  bs_kmers <- bootstrap_kmers(kmers, kmer_size)
  bs_class[[i]] <- classify_bs(bs_kmers, db)
}

consensus_bs_class(bs_class, db)

