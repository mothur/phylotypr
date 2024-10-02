#' Read in a FASTA-formatted file containing DNA sequences
#'
#' Given a
#' [standard FASTA-formatted file](https://en.wikipedia.org/wiki/FASTA_format),
#' `read_fasta` will read in the contents of the file and create a three column
#' data frame with columns for the sequence id, the sequence itself, and any
#' comments found in the header line for each sequence.
#'
#' @param file
#' Either a path to a file, a connection, or literal data (either a single
#' string or a raw vector) containing DNA sequences in the standard FASTA
#' format. There are no checks to determine whether the data are DNA or amino
#' acid sequences.
#'
#' Files ending in .gz, .bz2, .xz, or .zip will be automatically uncompressed.
#' Files starting with `http://`, `https://`, `ftp://`, or `ftps://` will be
#' automatically downloaded. Remote gz files can also be autom downloaded and
#' decompressed.
#'
#' @param degap
#' Logical value (default = TRUE) Removes gap characters from sequences
#' indicated by "." or "-"
#'
#' @note
#' The sequences in the FASTA file can have line breaks within them and
#' `read_fasta()` will put those separate lines into the same sequence
#'
#' @returns A data frame object with three columns. The `id` column will contain
#'   the non-space characters following the `>` in the header line of each
#'   sequence; the `sequence` column will contain the sequence; and the
#'   `comment` column will contain any text found after the first whitespace
#'   character on the header line.
#' @examples
#' temp <- tempfile()
#' write(">seqA\nATGCATGC\n>seqB\nTACGTACG", file = temp)
#' write(">seqC\nTCCGATGC", file = temp, append = TRUE)
#' write(">seqD B.ceresus UW85\nTCCGATGC", file = temp, append = TRUE)
#' write(">seq4\tE. coli K12\tBacteria;Proteobacteria;\nTCCGATGC",
#'   file = temp,
#'   append = TRUE
#' )
#' write(">seq_4\tSalmonella LT2\tBacteria;Proteobacteria;\nTCCGATGC",
#'   file = temp, append = TRUE
#' )
#' write(">seqE B.ceresus UW123\nTCCGATGC\nTCCGATGC",
#'   file = temp,
#'   append = TRUE
#' )
#'
#' sequence_df <- read_fasta(temp)
#'
#' @importFrom readr read_lines
#' @importFrom stringi stri_startswith_fixed stri_replace_first_regex stri_c
#' stri_replace_all_regex
#' @export
#'
read_fasta <- function(file, degap = TRUE) {
  fasta_data <- readr::read_lines(file)

  is_header <- stringi::stri_startswith_fixed(fasta_data, ">")
  header_lines <- fasta_data[is_header]

  id <- header_lines |>
    stringi::stri_replace_first_regex("^>(\\S*)\\s*.*", "$1")

  comment <- header_lines |>
    stringi::stri_replace_first_regex("^>\\S*\\s*(.*)", "$1")

  number <- cumsum(is_header)
  seq_lines <- fasta_data[!is_header]
  seq_number <- number[!is_header]

  sequence <- NULL

  if (degap) {
    sequence <- tapply(seq_lines, seq_number, stringi::stri_c, collapse = "") |>
      unname() |>
      stringi::stri_replace_all_regex("[.-]", "")
  } else {
    sequence <- tapply(seq_lines, seq_number, stringi::stri_c, collapse = "") |>
      unname()
  }

  data.frame(
    id = id,
    sequence = sequence,
    comment = comment
  )
}
