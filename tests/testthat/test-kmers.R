test_that("Can extract all possible 8-mers from a sequence", {
  x <- "ATGCGCTAGTAGCATGC"

  all_kmers <- get_all_kmers(x, kmer_size = 8)

  expected_kmers <- c(
    "ATGCGCTA", "TGCGCTAG", "GCGCTAGT", "CGCTAGTA",
    "GCTAGTAG", "CTAGTAGC", "TAGTAGCA", "AGTAGCAT",
    "GTAGCATG", "TAGCATGC"
  )

  all_kmers <- get_all_kmers(x)
  expect_equal(all_kmers, expected_kmers)
  expect_length(all_kmers, length(expected_kmers))


  all_kmers <- get_all_kmers(x, kmer_size = 9)

  expected_kmers <- c(
    "ATGCGCTAG", "TGCGCTAGT", "GCGCTAGTA", "CGCTAGTAG",
    "GCTAGTAGC", "CTAGTAGCA", "TAGTAGCAT", "AGTAGCATG",
    "GTAGCATGC"
  )

  expect_equal(all_kmers, expected_kmers)
  expect_length(all_kmers, length(expected_kmers))
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
  sequence <- "ATGCGCTAGTAGCATGC"

  kmers <- seq_to_base4(sequence) |> get_all_kmers()
  indices <- base4_to_index(kmers)

  detected <- detect_kmers(sequence)

  expect_equal(detected, indices)

  sequence <- "ATGCGCTAGTAGCATGCN"
  kmers <- seq_to_base4(sequence) |> get_all_kmers()
  indices <- base4_to_index(kmers)

  detected <- detect_kmers(sequence)

  expect_equal(detected, indices)


  sequence <- "ATGCGCTAGTAGCATGCN"
  kmers <- seq_to_base4(sequence) |> get_all_kmers(kmer_size = 7)
  indices <- base4_to_index(kmers)

  detected <- detect_kmers(sequence, kmer_size = 7)

  expect_equal(detected, indices)
})

test_that("Accurately detect kmers across multiple sequences", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC")
  base4_sequences <- seq_to_base4(sequences)

  expected <- vector(mode = "list", length = 2)
  expected[[1]] <- base4_to_index(get_all_kmers(base4_sequences[1], kmer_size))
  expected[[2]] <- base4_to_index(get_all_kmers(base4_sequences[2], kmer_size))

  detect_matrix <- detect_kmers_across_sequences(sequences, kmer_size)

  expect_equal(detect_matrix, expected)
})

test_that("Calcuate word specific priors", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  detect_list <- detect_kmers_across_sequences(sequences, kmer_size)

  # nolint start: commented_code_linter
  # 26 - all 3 = (3+0.5) / (1 + 3) =0.875
  # 29 - only 1 = 0.375
  # 30 - only 2 and 3 = 0.625
  # 64 - none = 0.125
  # nolint end: commented_code_linter

  priors <- calc_word_specific_priors(detect_list, kmer_size)

  expect_equal(priors[26], 0.875)
  expect_equal(priors[29], 0.375)
  expect_equal(priors[30], 0.625)
  expect_equal(priors[64], 0.125)
})


test_that("Calculate genus-specific conditional probabilities", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  genera <- c(1, 2, 2)

  detect_list <- detect_kmers_across_sequences(sequences, kmer_size)
  priors <- calc_word_specific_priors(detect_list, kmer_size)

  # nolint start: commented_code_linter
  # (m(wi) + Pi) / (M + 1)
  # 26 - all 3 = (c(1, 2)+0.5) / (c(1, 2) + 1) = 0.74 & 0.8333333
  # 29 - only 1 = (c(1, 0)+0.5) / (c(1, 2) + 1) = 0.7500000 0.1666667
  # 30 - only 2 and 3 = (c(0, 2)+0.5) / (c(1, 2) + 1) = 0.2500000 0.8333333
  # 64 - none = 0.125 = (c(0, 2)+0.5) / (c(1, 2) + 1) = 0.2500000 0.1666667
  # nolint end

  conditional_prob <- calc_genus_conditional_prob(
    detect_list,
    genera,
    priors
  )

  expect_equal(conditional_prob[26, ], log((c(1, 2) + 0.875) / (c(1, 2) + 1)))
  expect_equal(conditional_prob[29, ], log((c(1, 0) + 0.375) / (c(1, 2) + 1)))
  expect_equal(conditional_prob[30, ], log((c(0, 2) + 0.625) / (c(1, 2) + 1)))
  expect_equal(conditional_prob[64, ], log((c(0, 0) + 0.125) / (c(1, 2) + 1)))
})

test_that("Convert back and forth between genus names and indices", {
  genera_str <- c("A", "B", "B")
  genera_index <- c(1, 2, 2)

  expect_equal(genera_str_to_index(genera_str), genera_index)
  expect_equal(get_unique_genera(genera_str), c("A", "B"))
})

test_that("Create kmer database from sequences, taxonomy, and kmer size", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  genera <- c("A", "B", "B")

  db <- build_kmer_database(sequences, genera, kmer_size)

  expect_equal(
    db[["conditional_prob"]][26, ],
    log((c(1, 2) + 0.875) / (c(1, 2) + 1))
  )
  expect_equal(
    db[["conditional_prob"]][29, ],
    log((c(1, 0) + 0.375) / (c(1, 2) + 1))
  )
  expect_equal(
    db[["conditional_prob"]][30, ],
    log((c(0, 2) + 0.625) / (c(1, 2) + 1))
  )
  expect_equal(
    db[["conditional_prob"]][64, ],
    log((c(0, 0) + 0.125) / (c(1, 2) + 1))
  )

  expect_equal(db[["genera"]][1], "A")
  expect_equal(db[["genera"]][2], "B")
})

test_that("Bootstrap sample 1/kmer_size of kmers from unknowns", {
  kmers <- 1:100
  kmer_size <- 8
  expected_n_kmers <- as.integer(1 / 8 * 100)

  detected <- bootstrap_kmers(kmers, kmer_size)

  expect_length(detected, expected_n_kmers)
  expect_in(detected, kmers)
})


test_that("Classify boostrap sample works", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  genera <- c("A", "B", "B")

  db <- build_kmer_database(sequences, genera, kmer_size)
  unknown_kmers <- detect_kmers("ATGCGCTC", kmer_size)
  expected_classification <- 2

  detected_classification <- classify_bs(unknown_kmers, db)
  expect_equal(detected_classification, expected_classification)
})


test_that("Consensus classification of bootstrap subsamples", {
  db <- list()
  db[["genera"]] <- c("A;a;A", "A;a;B", "A;a;C", "A;b;A", "A;b;B", "A;b;C")

  bs_class <- c(1, 1, 1, 1, 4)

  expected <- list()
  expected[["taxonomy"]] <- c("A", "a", "A")
  expected[["confidence"]] <- c(1, 0.8, 0.8)

  observed <- consensus_bs_class(bs_class, db)

  expect_equal(observed, expected)
})


test_that("Return correct consensus taxonomy and confidence", {
  taxonomy <- c("A", "A", "A", "A", "A")
  expected <- list()
  expected[["frac"]] <- 1
  expected[["id"]] <- "A"
  observed <- get_consensus(taxonomy)
  expect_equal(observed, expected)


  taxonomy <- c("A;a", "A;a", "A;b", "A;b", "A;b")
  expected <- list()
  expected[["frac"]] <- 0.6
  expected[["id"]] <- "A;b"
  observed <- get_consensus(taxonomy)
  expect_equal(observed, expected)
})






test_that("Can classify a unknown sequence with a database", {
  kmer_size <- 3
  sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
  genera <- c("A", "B", "B")

  db <- build_kmer_database(sequences, genera, kmer_size)

  unknown_sequence <- "ATGCGCTC"
  expected <- list()
  expected[["taxonomy"]] <- "B"
  expected[["confidence"]] <- 1

  actual <- classify_sequence(
    unknown_sequence = unknown_sequence, database = db,
    kmer_size = kmer_size
  )
  expect_equal(actual, expected)
})
