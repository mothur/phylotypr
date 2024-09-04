library(microbenchmark)
set.seed(19760620)

fake_kmers <- sample(1:10, size = 50, replace = TRUE) + 20
n_kmers <- 64

k_table <- function() {
  kmer_table <- table(fake_kmers)
  kmer_values <- names(kmer_table) |> as.numeric()
  kmer_counts <- as.numeric(kmer_table)
  kmer_vector <- numeric(n_kmers)
  kmer_vector[kmer_values] <- kmer_counts
  kmer_vector
}

k_tabulate <- function() {
  tabulate(fake_kmers, nbins = n_kmers)
}

microbenchmark(k_table(), k_tabulate())
