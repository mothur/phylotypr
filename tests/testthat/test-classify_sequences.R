test_that("Can classify a unknown sequence with a database", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A", "B", "B")
  classification_database <- create_database(
    id = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")

  expected <- data.frame(
    id = as.character(1:3),
    level = c(1, 1, 1),
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
    level = 1,
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

  classification_database <- create_database(
    id = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1")

  expect_error(classify_sequences(
    sequences = sequences, database = db,
    kmer_size = kmer_size, ids = ids
  ))

  expected <- data.frame(
    id = as.character(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)),
    level = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
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

  classification_database <- create_database(
    id = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1", "seq_2", "seq_3")

  expected <- data.frame(
    id = as.character(c(
      "seq_1", "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    level = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
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

  classification_database <- create_database(
    id = 1:3,
    sequence = ref_sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

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
    level = c(1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4),
    taxonomy = c("A", "b", "C", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  expect_equal(actual, expected)
})


test_that("Can classify when a data frame is provided to data argument", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  classification_database <- create_database(
    id = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(s = sequences, i = ids)

  expected <- data.frame(
    id = as.character(c(
      "seq_1", "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    level = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
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


test_that("Test that strollur objects can be used in classify_sequences", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")
  reference <- data.frame(
    reference_names = "trainset9_pds",
    reference_vendor = "pds",
    reference_versions = "9",
    reference_usages = "classification",
    reference_notes = paste0("phylotypr v.", packageVersion("phylotypr")),
    reference_urls = paste0(
      "https://sourceforge.net/projects/rdp-classifier/",
      "files/RDP_Classifier_TrainingData/"
    ),
    reference_citation =
      "@article{doi:10.1128/AEM.00062-07,
        author = {Qiong Wang and George M. Garrity and James M. Tiedje and
          James R. Cole},
        title = {Naïve Bayesian Classifier for Rapid Assignment of
          rRNA Sequences into the New Bacterial Taxonomy},
        journal = {Applied and Environmental Microbiology},
        volume = {73},
        number = {16},
        pages = {5261-5267},
        year = {2007},
        doi = {10.1128/AEM.00062-07},

        URL = {https://journals.asm.org/doi/abs/10.1128/aem.00062-07},
        eprint = {https://journals.asm.org/doi/pdf/10.1128/aem.00062-07}
      }"
  )

  classification_database <- create_database(
    id = 1:3,
    sequence = sequences,
    taxonomy = genera,
    reference_data = reference
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(s = sequences, i = ids)
  strollur_object <- strollur::new_dataset("test")
  strollur::xdev_add_sequences(strollur_object, df,
    sequence_name = "i",
    sequence = "s", verbose = FALSE
  )

  strollur_object2 <- strollur::copy_dataset(strollur_object)

  expected <- strollur::xdev_assign_sequence_taxonomy_tidy(
    data = strollur_object2,
    table = data.frame(
      sequence_names = as.character(c(
        "seq_1", "seq_1", "seq_1", "seq_1",
        "seq_2", "seq_2", "seq_2", "seq_2",
        "seq_3", "seq_3", "seq_3", "seq_3"
      )),
      levels = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
      taxonomies = c(
        "A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"
      ),
      confidences = c(
        100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100
      )
    )
  )

  strollur::xdev_add_references(strollur_object2, table = reference)

  actual <- classify_sequences(
    data = strollur_object,
    database = db,
    kmer_size = kmer_size
  )

  expect_equal(
    strollur::xdev_report(
      strollur_object2,
      type = "references"
    )$reference_names,
    reference$reference_names
  )

  expect_equal(
    strollur::xdev_report(strollur_object, type = "references"),
    strollur::xdev_report(strollur_object2, type = "references")
  )

  # expect_equal_strollur_object2, strollur_object_ #need way to do this

  #####

  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(s = sequences, i = ids)

  strollur_object3 <- strollur::new_dataset("test")
  strollur::xdev_add_sequences(strollur_object3, df,
    sequence_name = "i",
    sequence = "s", verbose = FALSE
  )

  strollur_object4 <- strollur::copy_dataset(strollur_object3)

  expected <- strollur::xdev_assign_sequence_taxonomy_tidy(
    data = strollur_object4,
    table = data.frame(
      sequence_names = as.character(c(
        "seq_1", "seq_1", "seq_1", "seq_1",
        "seq_2", "seq_2", "seq_2", "seq_2",
        "seq_3", "seq_3", "seq_3", "seq_3"
      )),
      levels = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
      taxonomies = c(
        "A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"
      ),
      confidences = c(
        100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100
      )
    )
  )

  expected <- strollur::xdev_add_references(
    data = strollur_object4,
    table = reference
  )

  actual <- classify_sequences(
    data = strollur_object3,
    database = db,
    kmer_size = kmer_size
  )

  expect_equal(
    strollur::xdev_report(strollur_object4,
      type = "references"
    )$reference_names,
    "trainset9_pds"
  )

  expect_equal(
    strollur::xdev_report(strollur_object3, type = "references"),
    strollur::xdev_report(strollur_object4, type = "references")
  )

  # expect_equal_strollur_object, strollur_object2_ #need way to do this
})


test_that("Classify_sequences will error if strollur is not installed", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  classification_database <- create_database(
    id = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  ids <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(i = ids, s = sequences)
  strollur_object <- strollur::new_dataset("test")
  strollur::xdev_add_sequences(
    strollur_object, df,
    sequence_name = "i",
    sequence = "s", verbose = FALSE
  )

  local_mocked_bindings(require_namespace = function(...) FALSE)

  expect_error(classify_sequences(
    data = strollur_object,
    database = db,
    kmer_size = kmer_size
  ), "Cannot use this function without strollur installed.")
})
