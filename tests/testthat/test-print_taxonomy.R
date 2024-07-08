
test_that("Print out consesnsus taxonomy", {

  oscillospiraceae <- list(taxonomy = c("Bacteria", "Bacillota", "Clostridia",
                                        "Eubacteriales", "Oscillospiraceae"),
                           confidence = c(1.00, 1.00, 0.99, 0.99, 0.98))

  expected <- "Bacteria(100);Bacillota(100);Clostridia(99);Eubacteriales(99);Oscillospiraceae(98);Oscillospiraceae_unclassified(98)"

  tax_string <- print_taxonomy(oscillospiraceae, n_levels = 6)

  expect_equal(tax_string, expected)


  tax_string <- print_taxonomy(oscillospiraceae)
  expect_equal(tax_string, expected)


  bacteroidales <- list()
  bacteroidales[["taxonomy"]] <- c("Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales")
  bacteroidales[["confidence"]] <- c(1.00, 1.00, 0.97, 0.97)
  expected <- "Bacteria(100);Bacteroidota(100);Bacteroidia(97);Bacteroidales(97);Bacteroidales_unclassified(97);Bacteroidales_unclassified(97)"

  tax_string <- print_taxonomy(bacteroidales, n_levels = 6)

  expect_equal(tax_string, expected)

})
