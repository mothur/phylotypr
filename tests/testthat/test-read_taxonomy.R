test_that("reading in taxonomy data works", {
  temp <- tempfile()

  write("seqA\tA;B;C;", file = temp)
  write("seqB\tA;B; C;", file = temp, append = TRUE)
  write("seqC\tA; B;C;", file = temp, append = TRUE)
  write("seqD\tA;B;C", file = temp, append = TRUE)
  write("seqE\tA;B; C", file = temp, append = TRUE)
  write("seqF\tA; B;C", file = temp, append = TRUE)
  write("seq G\tA;B;C;", file = temp, append = TRUE)
  write("seq G1\tA;\"B\";C;", file = temp, append = TRUE)
  write("seq G2\tA;\"B\";\"C\";", file = temp, append = TRUE)
  write("seq G3\t\"A\";\"B\";\"C\";", file = temp, append = TRUE)
  write("seq H1\tA;\'B\';C;", file = temp, append = TRUE)
  write("seq H2\tA;\'B\';\'C\';", file = temp, append = TRUE)
  write("seq H3\t\'A\';\'B\';\'C\';", file = temp, append = TRUE)

  taxonomy_df <- read_taxonomy(temp)

  expected <- data.frame(
    id = c(
      "seqA", "seqB", "seqC", "seqD", "seqE", "seqF", "seq G",
      "seq G1", "seq G2", "seq G3", "seq H1", "seq H2", "seq H3"
    ),
    taxonomy = c("A;B;C")
  )

  expect_equal(taxonomy_df, expected)
})
