test_that("Apply filter to confidence score", {
  oscillospiraceae <- data.frame(
    sequence_name = "a",
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae",
      "Flintibacter"
    ),
    level = 1:6,
    confidence = c(100, 100, 99, 99, 98, 58)
  )

  filtered <- data.frame(
    sequence_name = "a",
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae"
    ),
    level = 1:5,
    confidence = c(100, 100, 99, 99, 98)
  )

  observed <- filter_taxonomy(oscillospiraceae)
  expect_equal(observed, filtered)

  observed <- filter_taxonomy(oscillospiraceae, min_confidence = 80)
  expect_equal(observed, filtered)

  # do sequence_name, taxonomy, and confidence columns exist
  filtered <- data.frame(
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae"
    ),
    level = 1:5,
    confidence = c(100, 100, 99, 99, 98)
  )
  expect_error(filter_taxonomy(filtered))

  # Is it a data.frame
  expect_error(filter_taxonomy(c()))
})

test_that("strollur objects can be filtered", {
  oscillospiraceae <- data.frame(
    sequence_name = "a",
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae",
      "Flintibacter"
    ),
    confidence = c(100, 100, 99, 99, 98, 58),
    level = 1:6
  )


  strollur_object <- strollur::new_dataset("data")

  strollur::xdev_assign_sequence_taxonomy_tidy(
    strollur_object,
    table = oscillospiraceae,
    sequence_name = "sequence_name",
    taxonomy = "taxonomy",
    level = "level",
    confidence = "confidence"
  )

  filtered <- data.frame(
    sequence_name = "a",
    taxonomy = c(
      "Bacteria", "Bacillota", "Clostridia",
      "Eubacteriales", "Oscillospiraceae"
    ),
    level = 1:5,
    confidence = c(100, 100, 99, 99, 98)
  )

  filtered_strollur_object <- strollur::new_dataset("data")

  strollur::xdev_assign_sequence_taxonomy_tidy(
    filtered_strollur_object,
    table = filtered,
    sequence_name = "sequence_name",
    taxonomy = "taxonomy",
    level = "level",
    confidence = "confidence"
  )

  filter_taxonomy(strollur_object)

  expect_equal(
    strollur::xdev_report(strollur_object, type = "sequence_taxonomy"),
    strollur::xdev_report(filtered_strollur_object, type = "sequence_taxonomy"))

  expect_true(strollur::is_equal(strollur_object, filtered_strollur_object))

})

