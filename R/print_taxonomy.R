#' Print taxonomy for an unknown sequence
#'
#' The `print_taxonomy()` will output the consensus taxonomy for an unknown
#' sequence with confidence scores for each taxonomic level and each taxonomic
#' level separated by semi-colons
#'
#' @inheritParams filter_taxonomy
#' @param n_levels   An integer indicating the number of taxonomic levels to
#'                   expect. If the number of observed levels is less than this
#'                   value, and `fill_levels` is `TRUE` then missing levels will
#'                   have "_unclassified" appended to the end of the last named
#'                   classification until the desired level is reached. If
#'                  `n_levels` is less than the observed number of levels, then
#'                   the string of the output will be truncated to `n_levels`.
#'                   If 0 is given then we will report the longest possible
#'                   taxonomy. Default is 0.
#' @param fill_levels A logical indicating whether or not to extend the
#'                    classifications of sequences that don't meet `n_levels`.
#'                    Default is `TRUE`.
#' @param show_confidence A logical indicating whether or not to include the
#'                        confidence scores. Default is `TRUE`
#'
#' @returns A vector of character strings indicating the classification at each
#'   taxonomic level with the corresponding confidence in parentheses (if
#'   `show_confidence` is `FALSE`. Each taxonomic level is separated by a
#'   semi-colon
#'
#' @examples
#' oscillospiraceae <- data.frame(
#'   id = c("seq_1", "seq_1", "seq_1", "seq_1", "seq_1"),
#'   taxonomy = c(
#'     "Bacteria", "Bacillota", "Clostridia",
#'     "Eubacteriales", "Oscillospiraceae"
#'   ),
#'   confidence = c(100, 100, 99, 99, 98)
#' )
#'
#' print_taxonomy(oscillospiraceae, n_levels = 6)
#' print_taxonomy(oscillospiraceae, n_levels = 4)
#' print_taxonomy(oscillospiraceae, n_levels = 6, show_confidence = FALSE)
#' print_taxonomy(oscillospiraceae, n_levels = 7, fill_levels = FALSE)

#' @export
print_taxonomy <- function(consensus, n_levels = 0, fill_levels = TRUE,
                           show_confidence = TRUE) {
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

  if (n_levels <= 0) {
    n_levels <- aggregate(taxonomy ~ id, consensus, length)$taxonomy |>
      max()
  }

  results <- if (show_confidence) {
    paste0(consensus$taxonomy, "(", consensus$confidence, ")")
  } else {
    consensus$taxonomy
  }

  results <- results |>
    aggregate(
      x = _, by = list(id = consensus$id),
      paste_taxonomy, n_levels = n_levels, fill_levels = fill_levels,
      show_confidence = show_confidence
    )

  colnames(results) <- c("id", "taxonomy")
  results
}


paste_taxonomy <- function(taxonomy, n_levels, fill_levels = TRUE,
                           show_confidence = TRUE) {
  observed_levels <- length(taxonomy)

  if ((n_levels > observed_levels) && fill_levels) {
    inserted_taxonomy <- if (show_confidence) {
      sub(
        pattern = "^(.*)(\\([0-9\\.]+\\))",
        replacement = "\\1_unclassified\\2",
        taxonomy[observed_levels]
      )
    } else {
      sub(
        pattern = "(.*)$",
        replacement = "\\1_unclassified",
        taxonomy[observed_levels]
      )
    }

    taxonomy[(observed_levels + 1):n_levels] <- inserted_taxonomy
  }

  taxonomy[1:n_levels] |>
    na.omit() |>
    paste0(collapse = ";")
}
