#' Build kmer database for classifying 16S rRNA and other gene sequences to
#' a genus when a kmer size is provided.
#'
#' @param sequences A vector of reference sequences for which we have
#'                  genus-level taxonomic information in the same order as the
#'                  value for genera.
#' @param genera    A vector of genus-level taxonomic information for reference
#'                  sequences in the same order as the value for sequences.
#'                  Ideally, taxonomic information will be provided back to the
#'                  domain level with each level separated by semicolons and no
#'                  spaces.
#' @param kmer_size The length of the nucleotide word to base our classification
#'                  on (default = 8)
#'
#' @return  A TBD object containing the genus level conditional probability of
#'          seeing each kmer in a given genus as well as the genus names
#' @export
build_kmer_database <- function(sequences, genera, kmer_size = 8) {

  genera_indices <- genera_str_to_index(genera)

  detected_kmers <- seq_to_base4(sequences) |>
    detect_kmers_across_sequences(kmer_size = kmer_size)

  priors <- calc_word_specific_priors(detected_kmers, kmer_size)

  cond_prob <- calc_genus_conditional_prob(detected_kmers, genera_indices,
                                           priors)
  genera_names <- get_unique_genera(genera)

  return(list(conditional_prob = cond_prob,
              genera = genera_names))
}


#' @noRd
get_all_kmers <- function(x, kmer_size = 8){

  n_kmers <- nchar(x) - kmer_size + 1

  seq_kmers <- character(length = n_kmers)
  for(i in seq_along(1:n_kmers)){
    seq_kmers[i] <- get_kmer(sequence = x, start = i, kmer_size = kmer_size)
  }
  return(seq_kmers)
}


#' @noRd
get_kmer <- function(sequence, start, kmer_size = 8) {

  if(start + kmer_size - 1 > nchar(sequence)) {
    stop("Cannot extract kmer beyond end of sequence")
  }

  substr(sequence, start, start + kmer_size - 1)

}


#' @noRd
seq_to_base4 <- function(sequence) {

  toupper(sequence) |>
  gsub(pattern = "[^ACGT]", replacement = "N", x = _)  |>
    chartr(old = "ACGT", new = "0123", x = _)

}


#' @noRd
base4_to_index <- function(base4_str) {
  # I want output to be indexed to start at position 1 rather than 0 so we're
  # adding 1 to all base10 values
  stats::na.omit(strtoi(base4_str, base = 4) + 1) |> as.numeric()
}


#' @noRd
detect_kmers <- function(sequence, kmer_size = 8) {

  kmers <- get_all_kmers(sequence, kmer_size)
  indices <- base4_to_index(kmers)
  return(indices)

}


#' @noRd
detect_kmers_across_sequences <- function(sequences, kmer_size = 8) {

  n_sequences <- length(sequences)
  kmer_list <- vector(mode = "list", length = n_sequences)

  for(i in seq_along(1:n_sequences)){
    kmer_list[[i]] <- detect_kmers(sequences[[i]], kmer_size = kmer_size)
  }

  return(kmer_list)
}


#' @noRd
calc_word_specific_priors <- function(detect_list, kmer_size) {

  priors <- detect_list |> unlist() |> tabulate(bin = _, nbins = 4^kmer_size)

  (priors + 0.5) / (length(detect_list) + 1)

}


#' @noRd
#(m(wi) + Pi) / (M + 1)
calc_genus_conditional_prob <- function(detect_list,
                                        genera, #needs to be an integer not char
                                        word_specific_priors) {

  genus_counts <- tabulate(genera)
  n_genera <- length(genus_counts)
  n_sequences <- length(genera)
  n_kmers <- length(word_specific_priors)

  genus_count <- matrix(0,
                        nrow = n_kmers,
                        ncol = n_genera)

  for(i in 1:n_sequences) {
    genus_count[detect_list[[i]], genera[i]] <- genus_count[detect_list[[i]], genera[i]] + 1
  }

  t(t(genus_count + word_specific_priors) / (genus_counts + 1))
}


#' @noRd
genera_str_to_index <- function(string) {

  factor(string) |> as.numeric()

}


#' @noRd
get_unique_genera <- function(string) {

  factor(string) |> levels()
}
