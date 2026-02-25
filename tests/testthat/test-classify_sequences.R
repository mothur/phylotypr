test_that("Can classify a unknown sequence with a database", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A", "B", "B")

  db <- build_kmer_database(sequences, genera, kmer_size)

  unknown_sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")

  expected <- data.frame(
    id = as.character(1:3),
    taxonomy = c("A", "B", "B"),
    confidence = c(100, 100, 100)
  )

  actual <- classify_sequences(
    unknown_sequences = unknown_sequences, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)


  unknown_sequences <- "AAAATCGAGTA" # test giving a vector of one sequence


  expected <- data.frame(
    id = c("1"),
    taxonomy = c("A"),
    confidence = c(100)
  )

  actual <- classify_sequences(
    unknown_sequences = unknown_sequences, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)

  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  db <- build_kmer_database(sequences, genera, kmer_size)

  unknown_sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")

  expected <- data.frame(
    id = as.character(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)),
    taxonomy = c("A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  actual <- classify_sequences(
    unknown_sequences = unknown_sequences, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)
})
