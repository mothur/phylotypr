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
