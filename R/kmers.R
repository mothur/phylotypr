get_all_kmers <- function(x, kmer_size = 8){

  n_kmers <- nchar(x) - kmer_size + 1
  sapply(1:n_kmers, get_kmer, sequence = x, kmer_size = kmer_size)

}


get_kmer <- function(sequence, start, kmer_size = 8) {

  if(start + kmer_size - 1 > nchar(sequence)) {
    stop("Cannot extract kmer beyond end of sequence")
  }

  substr(sequence, start, start + kmer_size - 1)

}

seq_to_base4 <- function(sequence) {

  toupper(sequence) |>
  gsub(pattern = "[^ACGT]", replacement = "N", x = _)  |>
    chartr(old = "ACGT", new = "0123", x = _)

}

base4_to_index <- function(base4_str) {
  # I want output to be indexed to start at position 1 rather than 0 so we're
  # adding 1 to all base10 values
  na.omit(strtoi(base4_str, base = 4) + 1) |> as.numeric()
}
