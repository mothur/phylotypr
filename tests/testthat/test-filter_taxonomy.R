test_that("Apply filter to confidence score", {
  oscillospiraceae <- list(
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae",
      "Flintibacter"
    ),
    confidence = c(100, 100, 99, 99, 98, 58)
  )

  filtered <- list(
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
})
