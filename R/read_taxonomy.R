#' Read in taxonomy files
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom stringi stri_replace_last_regex stri_replace_all_regex
#' @importFrom readr read_tsv cols col_character
read_taxonomy <- function(file) {

  taxonomy <- readr::read_tsv(file,
                  col_names = c("id", "taxonomy"),
                  col_types = readr::cols(.default = readr::col_character()))

  taxonomy$taxonomy <- stringi::stri_replace_last_regex(taxonomy$taxonomy, ";$", "")
  taxonomy$taxonomy <- stringi::stri_replace_all_regex(taxonomy$taxonomy, "; ", ";")
  as.data.frame(taxonomy)
}
