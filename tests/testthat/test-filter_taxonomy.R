test_that("Apply filter to confidence score", {
  oscillospiraceae <- data.frame(
    id = "a",
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae",
      "Flintibacter"
    ),
    confidence = c(100, 100, 99, 99, 98, 58)
  )

  filtered <- data.frame(
    id = "a",
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae"
    ),
    confidence = c(100, 100, 99, 99, 98)
  )

  observed <- filter_taxonomy(oscillospiraceae)
  expect_equal(observed, filtered)

  observed <- filter_taxonomy(oscillospiraceae, min_confidence = 80)
  expect_equal(observed, filtered)

  # do id, taxonomy, and confidence columns exist
  filtered <- data.frame(
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae"
    ),
    confidence = c(100, 100, 99, 99, 98)
  )
  expect_error(filter_taxonomy(filtered))

  # Is it a data.frame
  expect_error(filter_taxonomy(c()))
})
