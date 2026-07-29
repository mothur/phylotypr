library(phylotypr)
library(bench)

miseq <- read_fasta(phylotypr_example("miseq_sop.fasta.gz"))

db <- build_kmer_database(
  trainset9_pds$sequence,
  trainset9_pds$taxonomy
)

# kmers <- phylotypr:::detect_kmers(
#   sequence = miseq$sequence[1L]
# )
# kmers_int <- as.integer(kmers)

set.seed(19760620)
res <- mark(
  "R classify_bs" = phylotypr:::classify_bs(
    unknown_kmers = kmers,
    conditional_prob = db[[1L]]
  ),
  "C classify_bs" = phylotypr:::classify_c(
    unknown_kmers = kmers_int,
    conditional_prob = db[[1L]]
  ),
  "Rcpp matrix" = phylotypr:::classify_matrix(
    unknown_kmers = kmers_int,
    conditional_prob = db[[1L]]
  ),
  "Rcpp linear" = phylotypr:::classify_linear(
    unknown_kmers = kmers_int,
    conditional_prob = db[[1L]]
  ),
  iterations = 200L,
  check = TRUE
)

res

plot(res)
