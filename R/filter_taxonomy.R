#' Filter taxonomy
#'
#' The `filter_taxonomy()` function will filter a consensus taxonomy to remove
#' any taxonomic levels where the confidence score is below a `min_confidence`
#' level
#'
#' @inheritParams print_taxonomy
#' @param min_confidence A double value between 0 and 1 (default = 0.8). The
#'                       minimum fraction of bootstrap replicates that had
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
#' oscillospiraceae <- list(taxonomy = c("Bacteria", "Bacillota", "Clostridia",
#'                                       "Eubacteriales", "Oscillospiraceae",
#'                                       "Flintibacter"),
#'                          confidence = c(1.00, 1.00, 0.99, 0.99, 0.98, 0.58)
#'                        )
#'
#' filter_taxonomy(oscillospiraceae, min_confidence = 0.80)
#' @export

filter_taxonomy <- function(consensus, min_confidence = 0.80) {

  high_confidence <- which(consensus$confidence >= min_confidence)

  filtered <- list()
  filtered[["taxonomy"]] <- consensus[["taxonomy"]][high_confidence]
  filtered[["confidence"]] <- consensus[["confidence"]][high_confidence]

  return(filtered)
}
