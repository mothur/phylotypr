test_that("Can classify a unknown sequence with a database", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A", "B", "B")
  classification_database <- create_database(
    sequence_name = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")

  expected <- data.frame(
    sequence_name = as.character(1:3),
    level = c(1, 1, 1),
    taxonomy = c("A", "B", "B"),
    confidence = c(100, 100, 100)
  )

  # one processor
  actual <- classify_sequences(
    sequence = sequence, database = db, num_processors = 1,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)


  # two processors
  actual <- classify_sequences(
    sequence = sequence, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)


  sequence <- "AAAATCGAGTA" # test giving a vector of one sequence


  expected <- data.frame(
    sequence_name = c("1"),
    level = 1,
    taxonomy = c("A"),
    confidence = c(100)
  )

  actual <- classify_sequences(
    sequence = sequence, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)

  kmer_size <- 3
  sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  classification_database <- create_database(
    sequence_name = 1:3,
    sequence = sequence,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  sequence_name <- c("seq_1")

  expect_error(classify_sequences(
    sequence = sequence, database = db,
    kmer_size = kmer_size, sequence_name = sequence_name
  ))

  expected <- data.frame(
    sequence_name = as.character(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)),
    level = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
    taxonomy = c("A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  actual <- classify_sequences(
    sequence = sequence, database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)
})


test_that("Can classify a unknown sequence when names are provided", {
  kmer_size <- 3
  sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  classification_database <- create_database(
    sequence_name = 1:3,
    sequence = sequence,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  sequence_name <- c("seq_1", "seq_2", "seq_3")

  expected <- data.frame(
    sequence_name = as.character(c(
      "seq_1", "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    level = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
    taxonomy = c("A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  actual <- classify_sequences(
    sequence = sequence,
    sequence_name = sequence_name,
    database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)
})


test_that("classify_sequences with confidence filters the classification", {
  kmer_size <- 3
  ref_sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  classification_database <- create_database(
    sequence_name = 1:3,
    sequence = ref_sequence,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequence <- c("ATGATCGAGTC", "ATGCGCTCCCC", "ATGCGCTCCCC")
  sequence_name <- c("seq_1", "seq_2", "seq_3")

  set.seed(19760620)
  actual <- classify_sequences(
    sequence = sequence,
    sequence_name = sequence_name,
    database = db,
    kmer_size = kmer_size,
    min_confidence = 95
  )

  expected <- data.frame(
    sequence_name = as.character(c(
      "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    level = c(1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4),
    taxonomy = c("A", "b", "C", "A", "b", "C", "e", "A", "b", "C", "e"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  expect_equal(actual, expected)


  set.seed(19760620)
  actual <- classify_sequences(
    data.frame(sequence = sequence,
               sequence_name = sequence_name),
    database = db,
    kmer_size = kmer_size,
    min_confidence = 95
  )
  expect_equal(actual, expected)

})


test_that("Can classify when a data frame is provided to data argument", {
  kmer_size <- 3
  sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  classification_database <- create_database(
    sequence_name = 1:3,
    sequence = sequence,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  sequence_name <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(s = sequence, i = sequence_name)

  expected <- data.frame(
    sequence_name = as.character(c(
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
    sequence = "s",
    sequence_name = "i",
    database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)


  actual <- classify_sequences(
    data = df,
    sequence = s,
    sequence_name = i,
    database = db,
    kmer_size = kmer_size
  )

  expect_equal(actual, expected)


  expect_error(
    classify_sequences(
      data = df,
      sequence = a,
      sequence_name = i,
      database = db,
      kmer_size = kmer_size
    ),
    "Column named 'a' does not exist in 'df'"
  )

  expect_error(
    classify_sequences(
      data = df,
      sequence = s,
      sequence_name = j,
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
    name = "rdp_trainset9_pds",
    vendor = "rdp", # rdp, pds, silva, greengenes, ncbi
    version = "9_pds", # 10_pds, v8,
    usage = "classification",
    note = paste0("phylotypr v.", packageVersion("phylotypr")),
    documentation_url = paste0(
      "https://sourceforge.net/projects/rdp-classifier/",
      "files/RDP_Classifier_TrainingData/"
    ),
    citation =
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
    sequence_name = 1:3,
    sequence = sequences,
    taxonomy = genera,
    reference_data = reference
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  sequence_name <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(s = sequence, i = sequence_name)
  strollur_object <- strollur::new_dataset("test")
  strollur::xdev_add_sequences(strollur_object, df,
    sequence_name = "i",
    sequence = "s", verbose = FALSE
  )

  strollur_object2 <- strollur::copy_dataset(strollur_object)

  expected <- strollur::xdev_assign_sequence_taxonomy_tidy(
    data = strollur_object2,
    table = data.frame(
      sequence_name = as.character(c(
        "seq_1", "seq_1", "seq_1", "seq_1",
        "seq_2", "seq_2", "seq_2", "seq_2",
        "seq_3", "seq_3", "seq_3", "seq_3"
      )),
      level = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
      taxonomy = c(
        "A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"
      ),
      confidence = c(
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
      type = "resource_reference"
    )$reference_names,
    reference$reference_names
  )

  expect_equal(
    strollur::xdev_report(strollur_object, type = "resource_reference"),
    strollur::xdev_report(strollur_object2, type = "resource_reference")
  )

  expect_true(strollur::is_equal(strollur_object2, strollur_object))


  db <- build_kmer_database(classification_database, kmer_size)

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  sequence_name <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(s = sequences, i = sequence_name)

  strollur_object3 <- strollur::new_dataset("test")
  strollur::xdev_add_sequences(strollur_object3, df,
    sequence_name = "i",
    sequence = "s", verbose = FALSE
  )

  strollur_object4 <- strollur::copy_dataset(strollur_object3)

  expected <- strollur::xdev_assign_sequence_taxonomy_tidy(
    data = strollur_object4,
    table = data.frame(
      sequence_name = as.character(c(
        "seq_1", "seq_1", "seq_1", "seq_1",
        "seq_2", "seq_2", "seq_2", "seq_2",
        "seq_3", "seq_3", "seq_3", "seq_3"
      )),
      level = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
      taxonomy = c(
        "A", "b", "C", "d", "A", "b", "C", "e", "A", "b", "C", "e"
      ),
      confidence = c(
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
      type = "resource_reference"
    )$name,
    "rdp_trainset9_pds"
  )

  expect_equal(
    strollur::xdev_report(strollur_object3, type = "resource_reference"),
    strollur::xdev_report(strollur_object4, type = "resource_reference")
  )

  expect_true(strollur::is_equal(strollur_object3, strollur_object4))
})


test_that("Classify_sequences will error if strollur is not installed", {
  kmer_size <- 3
  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")

  classification_database <- create_database(
    sequence_name = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  db <- build_kmer_database(classification_database, kmer_size)

  sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  sequence_name <- c("seq_1", "seq_2", "seq_3")
  df <- data.frame(i = sequence_name, s = sequence)
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


test_that("Running multiple times with the same seed gives the same result and
          different result with different seeds", {

  db <- build_kmer_database(rdp_trainset9_pds)
  miseq <- read_fasta(phylotypr_example("miseq_sop.fasta.gz"))[c(1, 3, 5),]

  a <- miseq |>
    classify_sequences(database = db, seed = 1) |>
    print_taxonomy()

  a2 <- miseq |>
    classify_sequences(database = db, seed = 1) |>
    print_taxonomy()

  b <- miseq |>
    classify_sequences(database = db, seed = 2) |>
    print_taxonomy()

  expect_equal(a$taxonomy, a2$taxonomy)
  expect_false(all(a$taxonomy == b$taxonomy))

})
