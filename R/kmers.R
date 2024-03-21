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
