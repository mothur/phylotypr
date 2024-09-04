test_that("Read in fasta formatted data generates data frame", {
  temp <- tempfile()
  write(">seqA\nATGCATGC\n>seqB\nTACGTACG", file = temp)
  write(">seqC\nTCCGATGC", file = temp, append = TRUE)
  write(">seqD B.ceresus UW85\nTCCGATGC", file = temp, append = TRUE)
  write(">seq4\tE. coli K12\tBacteria;Proteobacteria;\nTCCGATGC",
    file = temp,
    append = TRUE
  )
  write(">seq_4\tSalmonella LT2\tBacteria;Proteobacteria;\nTCCGATGC",
    file = temp, append = TRUE
  )
  write(">seqE B.ceresus UW123\nTCCGATGC\nTCCGATGC",
    file = temp,
    append = TRUE
  )

  sequence_df <- read_fasta(temp)

  expected <- data.frame(
    id = c("seqA", "seqB", "seqC", "seqD", "seq4", "seq_4", "seqE"),
    sequence = c(
      "ATGCATGC", "TACGTACG", "TCCGATGC", "TCCGATGC", "TCCGATGC",
      "TCCGATGC", "TCCGATGCTCCGATGC"
    ),
    comment = c(
      "", "", "", "B.ceresus UW85",
      "E. coli K12\tBacteria;Proteobacteria;",
      "Salmonella LT2\tBacteria;Proteobacteria;", "B.ceresus UW123"
    )
  )

  expect_equal(sequence_df, expected)
})
