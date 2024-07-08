#' Print taxonomy for an unknown sequence
#'
#' The `print_taxonomy()` will output the consensus taxonomy for an unknown
#' sequence with confidence scores for each taxonomic level and each taxonomic
#' level separated by semi-colons
#'
#' @param consensus  A list object that contains two slots each with an equal
#'                   sized vector. The `taxonomy` vector contains the
#'                   classification at each taxonomic level and the `confidence`
#'                   vector contains the fraction of bootstraps that had the
#'                   specified classification
#' @param n_levels   An integer indicating the number of taxonomic levels to
#'                   expect. If the number of observed levels is less than this
#'                   value, then missing levels will have "_unclassified" to the
#'                   end of the last named classification
#'
#' @returns A character string indicating the classification at each taxonomic
#'   level with the corresponding confidence in parentheses. Each taxonomic
#'   level is separated by a semi-colon
#'
#' @examples
#' oscillospiraceae <- list(taxonomy = c("Bacteria", "Bacillota", "Clostridia",
#'                                       "Eubacteriales", "Oscillospiraceae"),
#'                           confidence = c(1.00, 1.00, 0.99, 0.99, 0.98)
#'                          )
#'
#' print_taxonomy(oscillospiraceae, n_levels = 6)

#' @export
print_taxonomy <- function(consensus, n_levels = 6) {

  original_levels <- length(consensus$taxonomy)
  given_levels <- original_levels

  while(given_levels < n_levels) {

    consensus$taxonomy[given_levels+1] <- paste(consensus$taxonomy[original_levels],
                                                "unclassified", sep = "_")
    consensus$confidence[given_levels+1] <- consensus$confidence[original_levels]
    given_levels <- given_levels + 1
  }

  pretty_confidence <- paste0("(", 100*consensus$confidence, ")")

  paste(consensus$taxonomy, pretty_confidence, sep = "", collapse = ";")

}

