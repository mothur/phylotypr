test_that("Print out consesnsus taxonomy", {
  classification <- data.frame(
    sequence_name = as.character(c(
      "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    level = c(1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5),
    taxonomy = c("A", "b", "C", "A", "b", "C", "e", "A", "b", "C", "e", "f"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )


  # We are truncating to a max of four and not filling anything less than four
  actual <- print_taxonomy(classification, n_levels = 3, fill_levels = FALSE)
  expected <- data.frame(
    sequence_name = c("seq_1", "seq_2", "seq_3"),
    taxonomy = c(
      "A(100);b(100);C(100)",
      "A(100);b(100);C(100)",
      "A(100);b(100);C(100)"
    )
  )
  expect_equal(actual, expected)

  actual <- print_taxonomy(classification, n_levels = 4, fill_levels = FALSE)
  expected <- data.frame(
    sequence_name = c("seq_1", "seq_2", "seq_3"),
    taxonomy = c(
      "A(100);b(100);C(100)",
      "A(100);b(100);C(100);e(100)",
      "A(100);b(100);C(100);e(100)"
    )
  )
  expect_equal(actual, expected)

  actual <- print_taxonomy(classification, n_levels = 0, fill_levels = FALSE)
  expected <- data.frame(
    sequence_name = c("seq_1", "seq_2", "seq_3"),
    taxonomy = c(
      "A(100);b(100);C(100)",
      "A(100);b(100);C(100);e(100)",
      "A(100);b(100);C(100);e(100);f(100)"
    )
  )
  expect_equal(actual, expected)

  actual <- print_taxonomy(classification, n_levels = 3, fill_levels = TRUE)
  expected <- data.frame(
    sequence_name = c("seq_1", "seq_2", "seq_3"),
    taxonomy = c(
      "A(100);b(100);C(100)",
      "A(100);b(100);C(100)",
      "A(100);b(100);C(100)"
    )
  )

  expect_equal(actual, expected)

  actual <- print_taxonomy(classification, n_levels = 4, fill_levels = TRUE)
  expected <- data.frame(
    sequence_name = c("seq_1", "seq_2", "seq_3"),
    taxonomy = c(
      "A(100);b(100);C(100);C_unclassified(100)",
      "A(100);b(100);C(100);e(100)",
      "A(100);b(100);C(100);e(100)"
    )
  )
  expect_equal(actual, expected)

  actual <- print_taxonomy(classification, n_levels = 0, fill_levels = TRUE)
  expected <- data.frame(
    sequence_name = c("seq_1", "seq_2", "seq_3"),
    taxonomy = c(
      "A(100);b(100);C(100);C_unclassified(100);C_unclassified(100)",
      "A(100);b(100);C(100);e(100);e_unclassified(100)",
      "A(100);b(100);C(100);e(100);f(100)"
    )
  )
  expect_equal(actual, expected)


  actual <- print_taxonomy(classification,
    n_levels = 0, fill_levels = TRUE,
    show_confidence = FALSE
  )
  expected <- data.frame(
    sequence_name = c("seq_1", "seq_2", "seq_3"),
    taxonomy = c(
      "A;b;C;C_unclassified;C_unclassified",
      "A;b;C;e;e_unclassified",
      "A;b;C;e;f"
    )
  )
  expect_equal(actual, expected)

  # do sequence_name, taxonomy, and confidence columns exist
  filtered <- data.frame(
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae"
    ),
    confidence = c(100, 100, 99, 99, 98)
  )
  expect_error(print_taxonomy(filtered))

  # Is it a data.frame
  expect_error(print_taxonomy(c()))
})


test_that("Paste consesnsus taxonomy under various situations", {
  taxonomy <- c("A(100)", "b(100)", "C(100)", "e(100)")
  observed <- paste_taxonomy(taxonomy, n_levels = 3, fill_levels = FALSE)
  expected <- "A(100);b(100);C(100)"
  expect_equal(observed, expected)

  observed <- paste_taxonomy(taxonomy, n_levels = 4, fill_levels = FALSE)
  expected <- "A(100);b(100);C(100);e(100)"
  expect_equal(observed, expected)

  observed <- paste_taxonomy(taxonomy, n_levels = 5, fill_levels = FALSE)
  expected <- "A(100);b(100);C(100);e(100)"
  expect_equal(observed, expected)


  observed <- paste_taxonomy(taxonomy, n_levels = 3, fill_levels = TRUE)
  expected <- "A(100);b(100);C(100)"
  expect_equal(observed, expected)

  observed <- paste_taxonomy(taxonomy, n_levels = 4, fill_levels = TRUE)
  expected <- "A(100);b(100);C(100);e(100)"
  expect_equal(observed, expected)

  observed <- paste_taxonomy(taxonomy, n_levels = 5, fill_levels = TRUE)
  expected <- "A(100);b(100);C(100);e(100);e_unclassified(100)"
  expect_equal(observed, expected)

  observed <- paste_taxonomy(taxonomy, n_levels = 5) # test default fill
  expected <- "A(100);b(100);C(100);e(100);e_unclassified(100)"
  expect_equal(observed, expected)

  taxonomy <- c("A(100)", "b(100)", "C(100)")
  observed <- paste_taxonomy(taxonomy, n_levels = 5, fill_levels = TRUE)
  expected <- "A(100);b(100);C(100);C_unclassified(100);C_unclassified(100)"


  taxonomy <- c("A(100)", "b(100)", "C(100)")
  observed <- paste_taxonomy(taxonomy, n_levels = 0, fill_levels = TRUE)
  expected <- "A(100);b(100);C(100)"
})

test_that("strollur objects can be printed", {
  classification <- data.frame(
    sequence_name = as.character(c(
      "seq_1", "seq_1", "seq_1",
      "seq_2", "seq_2", "seq_2", "seq_2",
      "seq_3", "seq_3", "seq_3", "seq_3", "seq_3"
    )),
    level = c(1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5),
    taxonomy = c("A", "b", "C", "A", "b", "C", "e", "A", "b", "C", "e", "f"),
    confidence = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  )

  strollur_object <- strollur::new_dataset("data")

  strollur::xdev_assign_sequence_taxonomy_tidy(
    strollur_object,
    table = classification,
    sequence_name = "sequence_name",
    taxonomy = "taxonomy",
    level = "level",
    confidence = "confidence"
  )

  actual <- print_taxonomy(strollur_object)

  expected <- data.frame(
    sequence_name = c("seq_1", "seq_2", "seq_3"),
    taxonomy = c(
      "A(100);b(100);C(100);C_unclassified(100);C_unclassified(100)",
      "A(100);b(100);C(100);e(100);e_unclassified(100)",
      "A(100);b(100);C(100);e(100);f(100)"
    )
  )
  expect_equal(actual, expected)



})
