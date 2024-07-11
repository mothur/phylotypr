## code to prepare `trainset19_df` dataset goes here

fasta <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.fasta"
taxonomy <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.tax"

fasta_df <- read_fasta(fasta)
genera <- read_taxonomy(taxonomy)

trainset19_df <- dplyr::inner_join(fasta_df, genera, by = "id")
trainset19_df <- trainset19[, c("id", "sequence", "taxonomy")]

usethis::use_data(trainset19_df, compress = "xz", overwrite = TRUE)


## code to prepare `trainset19_db` dataset goes here

trainset19_db <- build_kmer_database(trainset19_df$sequence,
                                     trainset19_df$taxonomy,
                                     kmer_size = 8)

usethis::use_data(trainset19_db, compress = "xz", overwrite = TRUE)
