#' Filter taxonomy
#'
#' The `filter_taxonomy()` function will filter a consensus taxonomy to remove
#' any taxonomic levels where the confidence score is below a `min_confidence`
#' level
#'
#' @param consensus  A data frame that contains at least three columns: the `id`
#'   column that contains an identifier for the sequence, the `taxonomy` column
#'   contains the classification at each taxonomic level, and the `confidence`
#'   column contains the percentage of bootstraps that had the specified
#'   classification
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
#'   id = c("seq_1", "seq_1", "seq_1", "seq_1", "seq_1", "seq_1"),
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
  # do id, taxonomy, and confidence columns exist
  if (!inherits(consensus, "data.frame")) {
    stop("`filter_taxonomy()` expects that consensus is a data.frame.")
  }
  if (!all(c("id", "taxonomy", "confidence") %in% colnames(consensus))) {
    message <- paste0(
      "`filter_taxonomy()` requires a data.frame with columns",
      " named id, taxonomy, and confidence."
    )
    stop(message)
  }
  consensus <- consensus[consensus$confidence >= min_confidence, ]
  rownames(consensus) <- NULL

  consensus
}
