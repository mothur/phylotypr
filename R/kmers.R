#' Build kmer database
#'
#' Build kmer database for classifying 16S rRNA and other gene sequences to
#' a genus when a kmer size is provided.
#'
#' @param sequences A vector of reference sequences for which we have
#'                  genus-level taxonomic information in the same order as the
#'                  value for genera.
#' @param genera    A character vector of genus-level taxonomic information for
#'                  reference sequences in the same order as the value for
#'                  sequences. Ideally, taxonomic information will be provided
#'                  back to the domain level with each level separated by
#'                  semicolons and no spaces.
#' @param kmer_size An integer indicating the length of the nucleotide word to
#'                  base our classification on (default = 8)
#' @param reference_data single row data frame object that has columns vendor,
#'   name, version, usage, notes, method_url, documentation_url, and citation.
#'   This is supplied with databases obtained from {phylotyprrefdata}. This
#'   information is only used if user supplied sequences are supplied as a
#'   {strollur} object.
#'
#' @return  A list object containing the genus level conditional probability
#'          (`conditional_prob`) of seeing each kmer in a given genus as well as
#'           the genus names (`genera`)
#'
#' @references
#' Wang Q, Garrity GM, Tiedje JM, Cole JR. Naive Bayesian classifier for rapid
#' assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ
#' Microbiol. 2007 Aug;73(16):5261-7.
#' doi:[10.1128/AEM.00062-07](https://pubmed.ncbi.nlm.nih.gov/17586664/)
#' PMID: 17586664; PMCID: PMC1950982.
#'
#' @examples
#' kmer_size <- 3
#' sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
#' genera <- c("A", "B", "B")
#'
#' build_kmer_database(sequences, genera, kmer_size)
#'
#' @export
build_kmer_database <- function(classification_database, kmer_size = 8) {
  if (!inherits(classification_database, "classification_database")) {
    stop("You supplied a database that is not of type classification_database")
  }

  genera_indices <- genera_str_to_index(
    classification_database$seq_taxonomy_data$taxonomy
  )

  detected_kmers <- DetectKmers(
    classification_database$seq_taxonomy_data$sequence,
    kmer_size = kmer_size
  )

  priors <- calc_word_specific_priors(detected_kmers, kmer_size)

  cond_prob <- calc_genus_conditional_prob(
    detected_kmers, genera_indices,
    priors
  )
  genera_names <- get_unique_genera(
    classification_database$seq_taxonomy_data$taxonomy
  )

  list(
    conditional_prob = cond_prob,
    genera = genera_names,
    reference_data = classification_database$reference_data
  )
}


#' Classify 16S rRNA gene sequence fragment
#'
#' The `classify_sequence()` function implements the Wang et al. naive Bayesian
#' classification algorithm for 16S rRNA gene sequences.
#'
#' @param sequence   A character object representing a DNA sequence that
#'                           needs to be classified.
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
#'                  replacement) from `sequence`. Higher values will
#'                  provide greater precision on the confidence score.
#'
#' @returns A list object of two vectors. One vector (`taxonomy`) is the
#'          taxonomic assignment for each level. The second vector
#'          (`confidence`) is the percentage of `num_bootstraps` that the
#'          classifier gave the same classification at that level.
#' @inherit build_kmer_database references
#' @export
#'
#' @examples
#' kmer_size <- 3
#' reference_sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
#' reference_genera <- c("A", "B", "B")
#'
#' db <- build_kmer_database(reference_sequences, reference_genera, kmer_size)
#' unknown_sequence <- "ATGCGCTC"
#'
#' classify_sequence(
#'   sequence = unknown_sequence,
#'   database = db,
#'   kmer_size = kmer_size
#' )
#'
classify_sequence <- function(sequence, database,
                              kmer_size = 8, num_bootstraps = 100) {
  .Deprecated("classify_sequences")
  classify_sequences(
    sequence = sequence, database = database,
    kmer_size = kmer_size, num_bootstraps = num_bootstraps
  )
}


#' @noRd
calc_word_specific_priors <- function(detect_list, kmer_size) {
  priors <- detect_list |>
    unlist() |>
    tabulate(bin = _, nbins = 4^kmer_size)

  (priors + 0.5) / (length(detect_list) + 1)
}


#' @noRd
# The forumula for calculating the conditional probability for each genus
# is equal to (m(wi) + Pi) / (M + 1) #nolint: commented_code_linter
# genera argument needs to be an integer not char
calc_genus_conditional_prob <- function(detect_list,
                                        genera,
                                        word_specific_priors) {
  genus_counts <- tabulate(genera)
  n_genera <- length(genus_counts)
  n_sequences <- length(detect_list)
  n_kmers <- length(word_specific_priors)

  kmer_genus_count <- matrix(0,
    nrow = n_kmers,
    ncol = n_genera
  )

  for (i in 1:n_sequences) {
    kmer_genus_count[detect_list[[i]], genera[i]] <-
      kmer_genus_count[detect_list[[i]], genera[i]] + 1
  }

  # this calls the Rcpp version of the function
  calculate_log_probability(
    kmer_genus_count,
    word_specific_priors,
    genus_counts
  ) |>
    t()
}


#' @noRd
genera_str_to_index <- function(string_matrix) {
  apply(string_matrix, 1, paste, collapse = ";") |>
    factor() |>
    as.numeric()
}


#' @noRd
get_unique_genera <- function(string_matrix) {
  apply(string_matrix, 1, paste, collapse = ";") |>
    factor() |>
    levels() |>
    strsplit(";") |>
    do.call(rbind, args = _) |>
    I()
}


#' @noRd
bootstrap_kmers <- function(kmers, kmer_size = 8) {
  n_kmers <- as.integer(length(kmers) / kmer_size)
  sample(kmers, n_kmers, replace = TRUE)
}


#' @noRd
#' @importFrom stringi stri_split_fixed stri_c
consensus_bs_class <- function(bs_class, genera) {
  lapply(
    bs_class,
    \(x) {
      # taxonomy <- genera[x]
      # # think about storing as already split in create_database...
      # stringi::stri_split_fixed(taxonomy, pattern = ";") |>
      #   get_consensus()

        get_consensus(genera[x, , drop = FALSE])

    }
  )
}

#' @noRd
get_consensus <- function(tax_table) {
  n_bootstraps <- nrow(tax_table)
  max_depth <- ncol(tax_table)

  max_count <- integer(length = max_depth)
  max_tax <- character(length = max_depth)

  prior_tax_count <- 0

  for (column in 1:max_depth) {

    tax_column <- tax_table[, column]

    column_factor <- factor(
      tax_column, levels = unique(tax_column)
    )

    tax_count <- tabulate(column_factor)
    max_index <- which.max(tax_count)
    max_count[column] <- tax_count[max_index]
    max_tax[column] <- levels(column_factor)[max_index]

    if(length(tax_count) == prior_tax_count){
      next
    }

    tax_table <- tax_table[tax_table[, column] == max_tax[column], ]
    prior_tax_count <- length(tax_count)
  }

  list(
    level = 1:max_depth,
    taxonomy = max_tax,
    confidence = 100 * max_count / n_bootstraps
  )
}
