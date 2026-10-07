#' Build custom database object for use with phylotypr
#'
#' The `create_database()` function builds a `classification_database` object
#' that can be used to generate the kmer database for downstream use in
#' `classify_sequences()`. This function is only necessary if you are building
#' your own database or modifying one of the existing databases provided in
#' `{phylotypr}` or `{phylotyprrefdata}`.
#'
#' @param sequence_name A character vector of identifiers for each sequence in
#'    the database. Should have the same number of values as found in the
#'    `sequence` and `taxonomy` vectors.
#' @param sequence A character vector of the sequences that are included in
#'    the database. Should have the same number of values as found in the
#'    `sequence_name` and `taxonomy` vectors.
#' @param taxonomy A character vector of taxonomy data for each sequence in
#'    the database. Should have the same number of values as found in the
#'    `sequence_name` and `sequence` vectors. The format of each value should be
#'    with each taxonomic level separated by a semicolon. Ideally each value in
#'    `taxonomy` will have the same number of levels.
#' @param reference_data single row data frame object that has columns vendor,
#'   name, version, usage, notes, method_url, documentation_url, and citation.
#'   This is supplied with databases obtained from \pkg{phylotyprrefdata}. This
#'   information is only used if user supplied sequences are supplied as a
#'   \pkg{strollur} object.
#'
#' @returns A `classification_database` object, which is a list containing
#'    the fields `seq_taxonomy_data` and `reference_data`. These are both
#'    data.frame objects.
#'
#' @examples
#' sequence_name <- c("seq_a", "seq_b", "seq_c")
#' sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
#' genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")
#'
#' create_database(
#'   sequence_name = sequence_name,
#'   sequence = sequence,
#'   taxonomy = genera
#' )
#'
#' @export
create_database <- function(sequence_name, sequence, taxonomy,
                            reference_data = NULL) {
  if (length(sequence_name) != length(sequence) ||
    length(sequence_name) != length(taxonomy)
  ) {
    stop("sequence_name, sequence, and taxonomy should all be the same length")
  }

  seq_taxonomy_data <- data.frame(sequence_name = sequence_name)
  seq_taxonomy_data["sequence"] <- clean_sequence(sequence)
  seq_taxonomy_data["taxonomy"] <- parse_taxonomy(taxonomy)

  structure(
    list(
      seq_taxonomy_data = seq_taxonomy_data,
      reference_data = reference_data
    ),
    class = "classification_database"
  )
}


parse_taxonomy <- function(taxonomy){

  split_taxonomy <- stringi::stri_split_fixed(taxonomy, pattern = ";")

  max_depth <- sapply(split_taxonomy, length) |> max()

  lapply(split_taxonomy, \(x) {
      # if number of levels < max_depth, add "" to get to max_depth
      if (length(x) < max_depth) {
        x <- c(x, rep(paste0(x[length(x)], "_unclassified"), max_depth - length(x)))
      } else {
        x
      }
    }) |>
    do.call(rbind, args = _) |>
    as.matrix() |>
    I()

}
