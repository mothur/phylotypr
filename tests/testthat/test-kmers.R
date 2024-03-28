test_that("Can extract all possible 8-mers from a sequence", {

  x <- "ATGCGCTAGTAGCATGC"

  all_kmers <- get_all_kmers(x, kmer_size = 8)

  expected_kmers <- c("ATGCGCTA", "TGCGCTAG", "GCGCTAGT", "CGCTAGTA",
                      "GCTAGTAG", "CTAGTAGC", "TAGTAGCA", "AGTAGCAT",
                      "GTAGCATG","TAGCATGC")

  all_kmers <- get_all_kmers(x)
  expect_equal(all_kmers, expected_kmers)
  expect_length(all_kmers, length(expected_kmers))


  all_kmers <- get_all_kmers(x, kmer_size = 9)

  expected_kmers <- c("ATGCGCTAG", "TGCGCTAGT", "GCGCTAGTA", "CGCTAGTAG",
                      "GCTAGTAGC", "CTAGTAGCA", "TAGTAGCAT", "AGTAGCATG",
                      "GTAGCATGC")

  expect_equal(all_kmers, expected_kmers)
  expect_length(all_kmers, length(expected_kmers))

})

test_that("Can extract specific kmer from a starting position and size", {

  x <- "ATGCGCTAGTAGCATGC"

  kmer <- get_kmer(x, 1, kmer_size = 8)
  expect_equal(kmer, "ATGCGCTA")

  kmer <- get_kmer(x, 5, kmer_size = 8)
  expect_equal(kmer, "GCTAGTAG")

  kmer <- get_kmer(x, 10, kmer_size = 8)
  expect_equal(kmer, "TAGCATGC")

  expect_error(get_kmer(x, 11, kmer_size = 8))

  kmer <- get_kmer(x, 5)
  expect_equal(kmer, "GCTAGTAG")


})


test_that("Conversion works between DNA sequence and quarternary", {

  x <- "ATGCGCTAGTAGCATGC"
  expected <- "03212130230210321"

  base4_seq <- seq_to_base4(x)
  expect_equal(base4_seq, expected)

  x <- "ATGCGCTRGTAGCATGC"
  expected <- "0321213N230210321"

  base4_seq <- seq_to_base4(x)
  expect_equal(base4_seq, expected)

  x <- tolower("ATGCGCTRGTAGCATGC")
  expected <- "0321213N230210321"

  base4_seq <- seq_to_base4(x)
  expect_equal(base4_seq, expected)

})

test_that("Generate base10 values from kmers in base4", {

  x <- "0000"
  expected <- 1
  actual <- base4_to_index(x)
  expect_equal(actual, expected)

  x <- "1000"
  expected <- 65
  actual <- base4_to_index(x)
  expect_equal(actual, expected)

  x <- "0123"
  expected <- 28
  actual <- base4_to_index(x)
  expect_equal(actual, expected)

  x <- c("0123", "1000", "0000")
  expected <- c(28, 65, 1)
  actual <- base4_to_index(x)
  expect_equal(actual, expected)

  x <- c("0123", "1000", "0000", "000N")
  expected <- c(28, 65, 1)
  actual <- base4_to_index(x)
  expect_equal(actual, expected)

})


test_that("Accurately detect kmers from a sequence", {

  sequence <- "03212130230210321"
  kmers <- get_all_kmers(sequence)
  indices <- base4_to_index(kmers)

  detected <- detect_kmers(sequence)

  expect_equal(length(detected[indices]), length(indices))
  expect_equal(length(detected[-indices]), 4^8 - length(indices))
  expect_equal(sum(detected[indices]), length(indices))


  sequence <- "03212130230210321N"
  kmers <- get_all_kmers(sequence)
  indices <- base4_to_index(kmers)

  detected <- detect_kmers(sequence)

  expect_equal(length(detected[indices]), length(indices))
  expect_equal(length(detected[-indices]), 4^8 - length(indices))
  expect_equal(sum(detected[indices]), length(indices))


  sequence <- "03212130230210321N"
  kmers <- get_all_kmers(sequence, kmer_size = 7)
  indices <- base4_to_index(kmers)

  detected <- detect_kmers(sequence, kmer_size = 7)

  expect_equal(length(detected[indices]), length(indices))
  expect_equal(length(detected[-indices]), 4^7 - length(indices))
  expect_equal(sum(detected[indices]), length(indices))


})

test_that("Accurately detect kmers across multiple sequences", {

  kmer_size <- 3
  sequences <- c("03212130", "03212131")

  expected <- matrix(0, nrow = 4^kmer_size, ncol = 2)
  expected[base4_to_index(get_all_kmers(sequences[1], kmer_size)), 1] <- 1
  expected[base4_to_index(get_all_kmers(sequences[2], kmer_size)), 2] <- 1

  detect_matrix <- detect_kmers_across_sequences(sequences, kmer_size)

  expect_equal(detect_matrix, expected)
})

test_that("Calcuate word specific priors", {

  kmer_size <- 3
  sequences <- c("03212130", "03212131", "03212131")
  detect_matrix <- detect_kmers_across_sequences(sequences, kmer_size)

  #26 - all 3 = (3+0.5) / (1 + 3) =0.875
  #29 - only 1 = 0.375
  #30 - only 2 and 3 = 0.625
  #64 - none = 0.125

  expected <- (apply(detect_matrix, 1, sum) + 0.5) / (1 + length(sequences))
  priors <- calc_word_specific_priors(detect_matrix)

  expect_equal(priors, expected)
  expect_equal(priors[26], 0.875)
  expect_equal(priors[29], 0.375)
  expect_equal(priors[30], 0.625)
  expect_equal(priors[64], 0.125)

})


test_that("Calculate genus-specific conditional probabilities", {

  kmer_size <- 3
  sequences <- c("03212130", "03212131", "03212131")
  genera <- c(1, 2, 2)

  detect_matrix <- detect_kmers_across_sequences(sequences, kmer_size)
  priors <- calc_word_specific_priors(detect_matrix)

  #(m(wi) + Pi) / (M + 1)
  #26 - all 3 = (c(1, 2)+0.5) / (c(1, 2) + 1) = 0.74 & 0.8333333
  #29 - only 1 = (c(1, 0)+0.5) / (c(1, 2) + 1) = 0.7500000 0.1666667
  #30 - only 2 and 3 = (c(0, 2)+0.5) / (c(1, 2) + 1) = 0.2500000 0.8333333
  #64 - none = 0.125 = (c(0, 2)+0.5) / (c(1, 2) + 1) = 0.2500000 0.1666667

  conditional_prob <- calc_genus_conditional_prob(detect_matrix,
                                                  genera,
                                                  priors)

  expect_equal(conditional_prob[26,], (c(1, 2)+0.875) / (c(1, 2) + 1))
  expect_equal(conditional_prob[29,], (c(1, 0)+0.375) / (c(1, 2) + 1))
  expect_equal(conditional_prob[30,], (c(0, 2)+0.625) / (c(1, 2) + 1))
  expect_equal(conditional_prob[64,], (c(0, 0)+0.125) / (c(1, 2) + 1))

})


test_that("Create kmer database from sequences, taxonomy, and kmer size", {

  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  genera <- c("A", "B", "B")

  db <- build_kmer_database(sequences, genera, kmer_size)

  expect_equal(db[["conditional_prob"]][26,], (c(1, 2)+0.875) / (c(1, 2) + 1))
  expect_equal(db[["conditional_prob"]][29,], (c(1, 0)+0.375) / (c(1, 2) + 1))
  expect_equal(db[["conditional_prob"]][30,], (c(0, 2)+0.625) / (c(1, 2) + 1))
  expect_equal(db[["conditional_prob"]][64,], (c(0, 0)+0.125) / (c(1, 2) + 1))

  expect_equal(db[["genera"]][1], "A")
  expect_equal(db[["genera"]][2], "B")

})

test_that("Convert back and forth between genus names and indices", {
  genera_str <- c("A", "B", "B")
  genera_index <- c(1, 2, 2)

  expect_equal(genera_str_to_index(genera_str), genera_index)
  expect_equal(get_unique_genera(genera_str), c("A", "B"))

})
