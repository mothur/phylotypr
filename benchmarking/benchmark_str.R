library(microbenchmark)
library(stringr)
library(stringi)

set.seed(19760620)

sequence <- sample(c("A", "T", "G", "C"), 1500, replace = T) |>
  paste(collapse = "")

# Calcualte the length of a string
# get_all_kmers / get_kmer
b_str_length <- function() {
  nchar(sequence)
}
i_str_length <- function() {
  stri_length(sequence)
}
r_str_length <- function() {
  str_length(sequence)
}

microbenchmark(b_str_length(), i_str_length(), r_str_length())


# Get substrings from a string
# get_kmer
n_kmers <- stringi::stri_length(sequence) - 8 + 1

b_substr <- function() {
  seq_kmers <- character(length = n_kmers)
  for (i in seq_along(1:n_kmers)) {
    seq_kmers[i] <- substr(sequence, i, i + 8 - 1)
  }
}

b_substring <- function() {
  substring(sequence, 1:n_kmers, 8:1500)
}
i_substring <- function() {
  stringi::stri_sub(sequence, 1:n_kmers, 8:1500)
}
r_substring <- function() {
  stringr::str_sub(sequence, 1:n_kmers, 8:1500)
}

microbenchmark(b_substr(), b_substring(), i_substring(), r_substrin)

# toupper and tolower
# seq_to_base4
sequence_upper <- sequence
sequence_lower <- tolower(sequence)

b_toupper <- function(x) {
  toupper(x)
}
b_tolower <- function(x) {
  tolower(x)
}
i_tolower <- function(x) {
  stringi::stri_trans_tolower(x)
}
i_toupper <- function(x) {
  stringi::stri_trans_toupper(x)
}
r_tolower <- function(x) {
  stringr::str_to_lower(x)
}
r_toupper <- function(x) {
  stringr::str_to_upper(x)
}

microbenchmark(
  b_toupper(sequence_upper),
  b_toupper(sequence_lower),
  b_tolower(sequence_upper),
  b_tolower(sequence_lower),
  i_toupper(sequence_upper),
  i_toupper(sequence_lower),
  i_tolower(sequence_upper),
  i_tolower(sequence_lower),
  r_toupper(sequence_upper),
  r_toupper(sequence_lower),
  r_tolower(sequence_upper),
  r_tolower(sequence_lower)
)


# want to replace non-atgc with a n character
sequence_good <- sequence
sequence_bad <- paste0(sequence_good, "R")

b_gsub <- function(x) {
  gsub(pattern = "[^ACGT]", replacement = "N", x = x)
}
i_replace <- function(x) {
  stringi::stri_replace_all_charclass(x, "[^ATGC]", "N")
}
r_replace <- function(x) {
  stringr::str_replace_all(x, "[^ATGC]", "N")
}

microbenchmark(
  b_gsub(sequence_good), b_gsub(sequence_bad),
  i_replace(sequence_good), i_replace(sequence_bad),
  r_replace(sequence_good), r_replace(sequence_bad)
)


# need to convert A C G T to 0 1 2 3

b_chartr <- function() {
  chartr("ACGT", "0123", sequence)
}
i_chartr <- function() {
  stringi::stri_trans_char(sequence, "ACGT", "0123")
}

microbenchmark(b_chartr(), i_chartr())
