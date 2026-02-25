#' Classify a vector of 16S rRNA gene sequence fragments
#'
#' The `classify_sequences()` function implements the Wang et al. naive Bayesian
#' classification algorithm for 16S rRNA gene sequences.
#'
#' @param unknown_sequences   A vector of character objects representing DNA
#'                            sequences that needs to be classified
#' @param database  A kmer database generated using [`build_kmer_database`]
#' @param kmer_size An integer value (default of 8) indicating the size of kmers
#'                  to use for classifying sequences. Higher values use more
#'                  RAM with potentially more specificity Lower values use
#'                  less RAM with potentially less specificity. Benchmarking
#'                  has found that the default of 8 provides the best
#'                  specificity with the lowest possible memory requirement and
#'                  fastest execution time.
#' @param num_bootstraps An integer value (default of 100). The value of
#'                  `num_bootstraps` is the number of randomizations to perform
#'                  where `1/kmer_size` of all kmers are sampled (without
#'                  replacement) from `unknown_sequence`. Higher values will
#'                  provide greater precision on the confidence score.
#'
#' @returns A list object of two vectors. One vector (`taxonomy`) is the
#'          taxonomic assignment for each level. The second vector
#'          (`confidence`) is the percentage of `num_bootstraps` that the
#'          classifier gave the same classification at that level
#' @inherit build_kmer_database references
#' @export
#'
#' @examples
#' kmer_size <- 3
#' sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
#' genera <- c("A", "B", "B")
#'
#' db <- build_kmer_database(sequences, genera, kmer_size)
#' unknown_sequence <- "ATGCGCTC"
#'
#' classify_sequence(
#'   unknown_sequence = unknown_sequence,
#'   database = db,
#'   kmer_size = kmer_size
#' )
classify_sequences <- function(unknown_sequences, database,
                               kmer_size = 8, num_bootstraps = 100) {
  lapply(seq_along(unknown_sequences), function(x) {
    list(
      id = as.character(x),
      classify_sequence(unknown_sequences[x],
        database = database,
        kmer_size = kmer_size, num_bootstraps = num_bootstraps
      )
    )
  }) |>
    lapply(X = _, as.data.frame) |>
    do.call(rbind, args = _)
}
