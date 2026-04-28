test_that("Database creation works", {
  actual <- create_database(
    id = rdp_trainset9_pds$seq_taxonomy_data$id,
    sequence = rdp_trainset9_pds$seq_taxonomy_data$sequence,
    taxonomy = rdp_trainset9_pds$seq_taxonomy_data$taxonomy,
    reference_data = rdp_trainset9_pds$reference_data
  )
  expected <- rdp_trainset9_pds

  expect_equal(actual, expected)
  expect_s3_class(actual, "classification_database")
})
