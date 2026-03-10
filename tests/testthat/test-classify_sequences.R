test_that("Can classify a unknown sequence with a database", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A", "B", "B")

  db <- build_kmer_database(sequences, genera, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")

  expected <- data.frame(
    id = as.character(1:3),
    taxonomy = c("A", "B", "B"),
    confidence = c(100, 100, 100)
  )

  actual <- classify_sequences(
    sequences = sequences, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)


  sequences <- "AAAATCGAGTA" # test giving a vector of one sequence


  expected <- data.frame(
    id = c("1"),
    taxonomy = c("A"),
    confidence = c(100)
  )

  actual <- classify_sequences(
    sequences = sequences, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)

  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  db <- build_kmer_database(sequences, genera, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1")

  expect_error(classify_sequences(
    sequences = sequences, database = db,
    kmer_size = kmer_size, ids = ids
  ))

  expected <- data.frame(
    id = as.character(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)),
    taxonomy = c("A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  actual <- classify_sequences(
    sequences = sequences, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)
})


test_that("Can classify a unknown sequence when names are provided", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  db <- build_kmer_database(sequences, genera, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1", "seq_2", "seq_3")

  expected <- data.frame(
    id = as.character(c(
      "seq_1", "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    taxonomy = c("A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  actual <- classify_sequences(
    sequences = sequences,
    ids = ids,
    database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)
})


test_that("classify_sequences with confidence filters the classification", {
  kmer_size <- 3
  ref_sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  db <- build_kmer_database(ref_sequences, genera, kmer_size)

  sequences <- c("ATGATCGAGTC", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1", "seq_2", "seq_3")

  set.seed(19760620)
  actual <- classify_sequences(
    sequences = sequences,
    ids = ids,
    database = db,
    kmer_size = kmer_size,
    min_confidence = 95
  )

  expected <- data.frame(
    id = as.character(c(
      "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    taxonomy = c("A", "b", "C", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  expect_equal(actual, expected)
})




test_that("Can classify when a data frame is provided to data argument", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  db <- build_kmer_database(sequences, genera, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(s = sequences, i = ids)

  expected <- data.frame(
    id = as.character(c(
      "seq_1", "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    taxonomy = c("A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  actual <- classify_sequences(
    data = df,
    sequences = "s",
    ids = "i",
    database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)


  actual <- classify_sequences(
    data = df,
    sequences = s,
    ids = i,
    database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)


  expect_error(
    classify_sequences(
      data = df,
      sequences = a,
      ids = i,
      database = db,
      kmer_size = kmer_size
    ),
    "Column named 'a' does not exist in 'df'"
  )

  expect_error(
    classify_sequences(
      data = df,
      sequences = s,
      ids = j,
      database = db,
      kmer_size = kmer_size
    ),
    "Column named 'j' does not exist in 'df'"
  )
})
