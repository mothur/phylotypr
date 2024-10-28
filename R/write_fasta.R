#' Write to a FASTA-formatted file
#'
#' Writes a data frame containing id, sequence, and comment columns,
#' `write_fasta` will write the data frame out to a
#' [standard FASTA-formatted file](https://en.wikipedia.org/wiki/FASTA_format).
#' The header will have a tab character between the sequence id and any
#' comments. There won't be a tab if there's no comment for the sequence. All
#' sequence data will be on a single line
#'
#' @param data_frame
#' A data frame object with three columns. The `id` column will contain
#' the non-space characters following the `>` in the header line of each
#' sequence; the `sequence` column will contain the sequence; and the
#' `comment` column will contain any text found after the first whitespace
#' character on the header line. The `comment` column is optional.
#'
#' @param file
#' Either a path to a file, a connection, or literal data (either a single
#' string or a raw vector) to write to a standard FASTA formatted file. There
#' are no checks to determine whether the data are DNA or amino acid sequences.
#'
#' Files ending in .gz, .bz2, .xz, or .zip will be automatically compressed.
#' Files starting with `http://`, `https://`, `ftp://`, or `ftps://` will be
#' automatically downloaded. Remote gz files can also be autom downloaded and
#' decompressed.
#'
#' If the value of `file` is `NULL` (default), the string will be written out
#' to the screen
#'
#' @returns
#' Sequence data is either written out to the screen (`file = NULL`) or to a
#' file.
#'
#' @examples
#' df_d <- data.frame(
#'   id = c("seqA", "seqB", "seqC"),
#'   sequence = c("ATGCATGC", "ATGCATGA", "ATGCATGT"),
#'   comment = c("comment 1", "", "comment 3")
#' )
#'
#' string_d <- write_fasta(df_d)
#' @export

write_fasta <- function(data_frame, file = NULL) {
  fasta <- paste0(">", data_frame$id, "\t", data_frame$comment, "\n",
    data_frame$sequence,
    collapse = "\n"
  ) |>
    gsub(pattern = "\t\n", replacement = "\n")

  if (!is.null(file)) {
    writeLines(fasta, file)
  } else {
    fasta
  }
}
