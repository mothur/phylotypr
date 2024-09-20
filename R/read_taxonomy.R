#' Read in taxonomy files
#'
#' Read a
#' [mothur-formatted taxonomy file](https://mothur.org/wiki/taxonomy_file/)
#' into R as a data frame
#'
#' @param file
#' Either a path to a file, a connection, or literal data
#' (either a single string or a raw vector) containing the
#' sequence id and the taxonomy information for each sequence.
#'
#' Files ending in .gz, .bz2, .xz, or .zip will be automatically uncompressed.
#' Files starting with ⁠http://⁠, ⁠https://⁠, ⁠ftp://⁠, or ⁠ftps://⁠ will be
#' automatically downloaded. Remote gz files can also be autom downloaded and
#' decompressed.
#'
#' @returns
#' A data frame with two columns. The `id` column contains a name for each
#' sequence and the `taxonomy` column, which contains the taxonomy for each
#' sequence. The string in the `taxonomy` column is a series of taxonomic names
#' separated by semi-colons. The string does not have a semi-colon at the end of
#' the sequence
#'
#' @note
#' There are no checks to insure that each sequence has a unique id value. It is
#' also assumed that each sequence has the same number of taxonomic levels
#' represented in the second column of the input file.
#'
#' @examples
#' temp <- tempfile()
#' write("seqA\tA;B;C;", file = temp)
#' write("seqB\tA;B; C;", file = temp, append = TRUE)
#' write("seqC\tA; B;C;", file = temp, append = TRUE)
#' write("seqD\tA;B;C", file = temp, append = TRUE)
#' write("seqE\tA;B; C", file = temp, append = TRUE)
#' write("seqF\tA; B;C", file = temp, append = TRUE)
#' write("seq G\tA;B;C;", file = temp, append = TRUE)
#'
#' read_taxonomy(temp)
#' @importFrom stringi stri_replace_last_regex stri_replace_all_regex
#' @importFrom readr read_tsv cols col_character
#' @export
read_taxonomy <- function(file) {
  taxonomy <- readr::read_tsv(file,
    col_names = c("id", "taxonomy"),
    col_types = readr::cols(.default = readr::col_character())
  )

  taxonomy$taxonomy <- stringi::stri_replace_last_regex(
    taxonomy$taxonomy,
    ";$",
    ""
  )

  taxonomy$taxonomy <- stringi::stri_replace_all_regex(
    taxonomy$taxonomy,
    "; ",
    ";"
  )

  taxonomy$taxonomy <- stringi::stri_replace_all_regex(
    taxonomy$taxonomy,
    "[\"\']",
    ""
  )

  as.data.frame(taxonomy)
}
