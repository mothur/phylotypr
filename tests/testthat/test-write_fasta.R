test_that("Output fasta-formatted data from a sequence data frame", {
  df_a <- data.frame(
    id = "seqA",
    sequence = "ATGCATGC",
    comment = ""
  )

  string_a <- write_fasta(df_a)
  expect_equal(string_a, ">seqA\nATGCATGC")


  df_b <- data.frame(
    id = "seqA",
    sequence = "ATGCATGC",
    comment = "this is a comment"
  )

  string_b <- write_fasta(df_b)
  expect_equal(string_b, ">seqA\tthis is a comment\nATGCATGC")


  df_c <- data.frame(
    id = "seqA",
    sequence = "ATGCATGC"
  )

  string_c <- write_fasta(df_c)
  expect_equal(string_c, ">seqA\nATGCATGC")


  df_d <- data.frame(
    id = c("seqA", "seqB", "seqC"),
    sequence = c("ATGCATGC", "ATGCATGA", "ATGCATGT"),
    comment = c("comment 1", "", "comment 3")
  )

  string_d <- write_fasta(df_d, file = NULL)
  expect_equal(
    string_d,
    ">seqA\tcomment 1\nATGCATGC\n>seqB\nATGCATGA\n>seqC\tcomment 3\nATGCATGT"
  )

  temp <- tempfile(fileext = ".fasta")
  write_fasta(df_d, file = temp)
  file_d <- read_fasta(temp)
  expect_equal(file_d, df_d)
})
