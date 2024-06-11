#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom readr read_lines
#' @importFrom stringi stri_startswith_fixed stri_replace_first_regex
read_fasta <- function(file) {

  fasta_data <- readr::read_lines(file)

  is_header <- stringi::stri_startswith_fixed(fasta_data, ">")
  header_lines <- fasta_data[is_header]

  id <- header_lines |>
    stringi::stri_replace_first_regex("^>(\\w*)\\s*.*", "$1")

  comment <- header_lines |>
    stringi::stri_replace_first_regex("^>\\w*\\s*(.*)", "$1")

  number <- cumsum(is_header)
  seq_lines <- fasta_data[!is_header]
  seq_number <- number[!is_header]

  sequence <- tapply(seq_lines, seq_number, stringi::stri_c, collapse = "") |>
    unname()

  data.frame(id = id,
             sequence = sequence,
             comment = comment)
}
