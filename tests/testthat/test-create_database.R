test_that("Database creation works", {
  actual <- create_database(
    sequence_name = rdp_trainset9_pds$seq_taxonomy_data$sequence_name,
    sequence = rdp_trainset9_pds$seq_taxonomy_data$sequence,
    taxonomy = apply(rdp_trainset9_pds$seq_taxonomy_data$taxonomy, 1, stri_c, collapse = ";"),
    reference_data = rdp_trainset9_pds$reference_data
  )
  expected <- rdp_trainset9_pds

  expect_equal(actual, expected)
  expect_s3_class(actual, "classification_database")
})

test_that("Database creation fails when fewer taxonomies than sequences", {

  expect_error(
    create_database(
      sequence_name = rdp_trainset9_pds$seq_taxonomy_data$sequence_name,
      sequence = rdp_trainset9_pds$seq_taxonomy_data$sequence,
      taxonomy = rdp_trainset9_pds$seq_taxonomy_data$taxonomy[-1],
      reference_data = rdp_trainset9_pds$reference_data
    ),
    "sequence_name, sequence, and taxonomy should all be the same length"
  )
})


test_that("We appropriately separate taxonomic levels by semicolon and add
  unclassified label", {

    taxonomy <- c("A;b;C;d", "A;b;C", "A;b;C", "A;b;C", "A;e;C")

    expected <- I(matrix(c("A", "A", "A", "A", "A",
                         "b", "b", "b", "b", "e",
                         "C", "C", "C", "C", "C",
                         "d", "C_unclassified", "C_unclassified",
                          "C_unclassified", "C_unclassified"), nrow = 5))


    observed <- parse_taxonomy(taxonomy)
    expect_equal(observed, expected)

  })

test_that("If only one taxonomic level, keep as matrix", {

  sequences <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
  genera <- c("A", "B", "B")
  classification_database <- create_database(
    sequence_name = 1:3,
    sequence = sequences,
    taxonomy = genera
  )
  seq_taxonomy_data <- data.frame(sequence_name = 1:3)
  seq_taxonomy_data["sequence"] <- sequences
  seq_taxonomy_data["taxonomy"] <- I(matrix(genera, ncol = 1))

  expected <- list(seq_taxonomy_data = seq_taxonomy_data,
                   reference_data = NULL)
  class(expected) <- "classification_database"

  expect_equal(classification_database$seq_taxonomy_data$sequence_name, expected$seq_taxonomy_data$sequence_name)
  expect_equal(classification_database$seq_taxonomy_data$sequence, expected$seq_taxonomy_data$sequence)
  expect_equal(classification_database$seq_taxonomy_data$taxonomy, expected$seq_taxonomy_data$taxonomy)
  expect_equal(classification_database, expected)
})
