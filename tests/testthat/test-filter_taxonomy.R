test_that("Apply filter to confidence score", {

  oscillospiraceae <- list(taxonomy = c("Bacteria", "Bacillota", "Clostridia",
                                        "Eubacteriales", "Oscillospiraceae",
                                        "Flintibacter"),
                           confidence = c(1.00, 1.00, 0.99, 0.99, 0.98, 0.58))

  filtered <- list(taxonomy = c("Bacteria", "Bacillota", "Clostridia",
                                "Eubacteriales", "Oscillospiraceae"),
                   confidence = c(1.00, 1.00, 0.99, 0.99, 0.98))

  observed <- filter_taxonomy(oscillospiraceae, min_confidence = 0.80)
  expect_equal(observed, filtered)

})
