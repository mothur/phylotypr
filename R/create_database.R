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
#' @inheritParams build_kmer_database reference_data
#'
#' @returns A `classification_database` object, which is a list containing
#'    the fields `seq_taxonomy_data` and `reference_data`. These are both
#'    data.frame objects.
#'
#' @examples
#'   sequence_name <- c("seq_a", "seq_b", "seq_c")
#'   sequence <- c("AAAATCGAGTA", "ATGCGCTCCCC", "ATGCGCTCCCC")
#'   genera <- c("A;b;C;d", "A;b;C;e", "A;b;C;e")
#'
#'   create_database(
#'     sequence_name = sequence_name,
#'     sequence = sequence,
#'     taxonomy = genera)
#'
#' @export
create_database <- function(sequence_name, sequence, taxonomy,
                            reference_data = NULL) {

  if(length(sequence_name) != length(sequence) ||
     length(sequence_name) != length(taxonomy)) {
    stop("sequence_name, sequence, and taxonomy should all be the same length")
  }

  structure(
    list(
      seq_taxonomy_data = data.frame(
        sequence_name = sequence_name,
        sequence = sequence,
        taxonomy = taxonomy
      ),
      reference_data = reference_data
    ),
    class = "classification_database"
  )
}
