#' Filter taxonomy
#'
#' The `filter_taxonomy()` function will filter a consensus taxonomy to remove
#' any taxonomic levels where the confidence score is below a `min_confidence`
#' level
#'
#' @param consensus  A data frame that contains at least three columns: the
#'   `sequence_name` column that contains an identifier for the sequence, the
#'   `taxonomy` column contains the classification at each taxonomic level, and
#'   the `confidence` column contains the percentage of bootstraps that had the
#'   specified classification. The `level` column is optional and contains the
#'   numeric level of the level (e.g., Kingdom = 1, Genus = 6)
#' @param min_confidence A double value between 0 and 100 (default = 80). The
#'   minimum percentage of bootstrap replicates that had the same
#'   classification. Any confidence score below this value will have the
#'   corresponding taxonomy removed
#'
#' @return A data frame object that is filtered to remove low confidence
#'   taxonomies. One column, `taxonomy`, contains the taxonomy at each taxonomic
#'   level and the other column, `confidence` contains the confidence score for
#'   that taxonomy. There will be no taxonomies or confidence scores below
#'   `min_confidence`
#'
#' @examples
#' oscillospiraceae <- data.frame(
#'   sequence_name = c("seq_1", "seq_1", "seq_1", "seq_1", "seq_1", "seq_1"),
#'   level = 1:6,
#'   taxonomy = c(
#'     "Bacteria", "Bacillota", "Clostridia",
#'     "Eubacteriales", "Oscillospiraceae",
#'     "Flintibacter"
#'   ),
#'   confidence = c(100, 100, 99, 99, 98, 58)
#' )
#'
#' filter_taxonomy(oscillospiraceae, min_confidence = 80)
#' @export
filter_taxonomy <- function(consensus, min_confidence = 80) {
  if (!inherits(consensus, c("data.frame", "strollur"))) {
    stop("`filter_taxonomy()` expects that consensus is a data.frame or strollur object.")
  }
  UseMethod("filter_taxonomy", consensus)
  }

#' @export
filter_taxonomy.data.frame <- function(consensus, min_confidence = 80) {
  do_filter(consensus, min_confidence)

}

#' @export
filter_taxonomy.strollur <- function(consensus, min_confidence = 80) {

  df <- strollur::xdev_report(consensus, "sequence_taxonomy")
  filtered_df <- do_filter(df, min_confidence)

  strollur::xdev_assign_sequence_taxonomy_tidy(
    consensus,
    table = filtered_df,
    sequence_name = "sequence_name",
    taxonomy = "taxonomy",
    level = "level",
    confidence = "confidence"
  )

}


do_filter <- function(consensus, min_confidence) {
  if (!all(c("sequence_name", "taxonomy", "level", "confidence") %in% colnames(consensus))) {
    message <- paste0(
      "`filter_taxonomy()` requires a data.frame with columns",
      " named sequence_name, taxonomy, level, and confidence."
    )
    stop(message)
  }
  consensus <- consensus[consensus$confidence >= min_confidence, ]
  rownames(consensus) <- NULL

  consensus
  }
