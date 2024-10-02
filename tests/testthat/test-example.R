test_that("NULL path returns names of example files", {
  example_files <- phylotypr_example()
  expect_true("miseq_sop.fasta.gz" %in% example_files)
  expect_true(is.character(example_files))
})

test_that("providing example file name returns full path", {
  expect_true(file.exists(phylotypr_example("miseq_sop.fasta.gz")))
})
