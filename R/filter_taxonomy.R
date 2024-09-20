#' Filter taxonomy
#'
#' The `filter_taxonomy()` function will filter a consensus taxonomy to remove
#' any taxonomic levels where the confidence score is below a `min_confidence`
#' level
#'
#' @inheritParams print_taxonomy
#' @param min_confidence A double value between 0 and 100 (default = 80). The
#'                       minimum percentage of bootstrap replicates that had
#'                       the same classification. Any confidence score below
#'                       this value will have the corresponding taxonomy
#'                       removed
#'
#' @return A list object containing two equally sized vectors that are filtered
#'         to remove low confidence taxonomies. One vector, `taxonomy`,
#'         contains the taxonomy at each taxonomic level and the other vector,
#'         `confidence` contains the confidence score for that taxonomy. There
#'         will be no taxonomies or confidence scores below `min_confidence`
#'
#' @examples
#' oscillospiraceae <- list(
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
  high_confidence <- which(consensus$confidence >= min_confidence)

  filtered <- list()
  filtered[["taxonomy"]] <- consensus[["taxonomy"]][high_confidence]
  filtered[["confidence"]] <- consensus[["confidence"]][high_confidence]

  return(filtered)
}
