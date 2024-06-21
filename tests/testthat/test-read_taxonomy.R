test_that("reading in taxonomy data works", {

  temp <- tempfile()

  write("seqA\tA;B;C;", file = temp)
  write("seqB\tA;B; C;", file = temp, append = TRUE)
  write("seqC\tA; B;C;", file = temp, append = TRUE)
  write("seqD\tA;B;C", file = temp, append = TRUE)
  write("seqE\tA;B; C", file = temp, append = TRUE)
  write("seqF\tA; B;C", file = temp, append = TRUE)
  write("seq G\tA;B;C;", file = temp, append = TRUE)

  taxonomy_df <- read_taxonomy(temp)

  expected <- data.frame(
    id = c("seqA", "seqB", "seqC", "seqD", "seqE", "seqF", "seq G"),
    taxonomy = c("A;B;C")

  )

  expect_equal(taxonomy_df, expected)

})
