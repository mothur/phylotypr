#library(dplyr)

bacteroidales <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTTAAGTCAGTGGTCAAATTGAGGGGCTCAACCCCTTCCCGCCATTGAAACTGGCGATCTTGAGTGGAAGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATGCCGGCTTCCTACTGACGCTCATGCACGAAAGTGTGGGTAACGAACAGG"

num_bootstraps <- 100
kmer_size <- 8


db <- build_kmer_database(trainset19_df$sequence,
                          trainset19_df$taxonomy,
                          kmer_size = 8)

consensus <- classify_sequence(unknown = bacteroidales, database = db,
                               num_bootstraps = num_bootstraps,
                               kmer_size = kmer_size)

filtered <- filter_taxonomy(consensus)
print_taxonomy(filtered)

