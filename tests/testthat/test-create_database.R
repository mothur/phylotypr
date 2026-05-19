test_that("Database creation works", {
  actual <- create_database(
    sequence_name = rdp_trainset9_pds$seq_taxonomy_data$sequence_name,
    sequence = rdp_trainset9_pds$seq_taxonomy_data$sequence,
    taxonomy = rdp_trainset9_pds$seq_taxonomy_data$taxonomy,
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
