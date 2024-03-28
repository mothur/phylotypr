library(tidyverse)

genera <- read_tsv("trainset19_072023.rdp/trainset19_072023.rdp.tax",
                   col_names = c("accession", "taxonomy"))

fasta_data <- scan("trainset19_072023.rdp/trainset19_072023.rdp.fasta",
                   sep = "\n",
                   what = character(), quiet = TRUE)

sequence_names <- fasta_data[seq(1, length(fasta_data), by = 2)] |>
  str_replace_all(pattern = "^>([^\t]*)\t.*", "\\1")

sequences <- fasta_data[seq(2, length(fasta_data), by = 2)]

seq_table <- tibble(accession = sequence_names, sequence = sequences) |>
  inner_join(genera, by = "accession")

build_kmer_database(seq_table$sequence, seq_table$taxonomy, kmer_size = 8)
