#' Classify a vector of 16S rRNA gene sequence fragments
#'
#' The `classify_sequences()` function implements the Wang et al. naive Bayesian
#' classification algorithm for 16S rRNA gene sequence.
#'
#' @param data      A data frame containing a column of `sequence_name` and
#'   `sequence`. If provided column names need to be provided to `sequence_name`
#'   and `sequence`.
#' @param sequence  A vector of character objects representing DNA
#'                  sequence that need to be classified. If `data` is provided
#'                  then only the name of the column containing sequence needs
#'                  to be provided, quotation marks are optional. This argument
#'                  can be ignored when using a `{strollur}` object
#' @param sequence_name       A vector of character objects containing the
#'                  sequence name that need to be classified. Optional. Defaults
#'                  to the position of the sequence in the vector of sequence.
#'                  If `data` is provided then only the name of the column
#'                  containing the sequence name needs to be provided, quotation
#'                  marks are optional. This argument can be ignored when using
#'                  a `{strollur}` object
#' @param database  A kmer database generated using [`build_kmer_database`]
#' @param kmer_size An integer value (default of 8) indicating the size of kmers
#'                  to use for classifying sequence. Higher values use more
#'                  RAM with potentially more specificity Lower values use
#'                  less RAM with potentially less specificity. Benchmarking
#'                  has found that the default of 8 provides the best
#'                  specificity with the lowest possible memory requirement and
#'                  fastest execution time.
#' @param num_bootstraps An integer value (default of 100). The value of
#'                  `num_bootstraps` is the number of randomizations to perform
#'                  where `1/kmer_size` of all kmers are sampled (without
#'                  replacement) from `unknown_sequence`. Higher values will
#'                  provide greater precision on the confidence score.
#' @param min_confidence An integer value (default of 80). The value of
#'                  `min_confidence` represents the minimum number of bootstraps
#'                  that must have a classification that agree with each other.
#' @param seed An integer value (default of 1). The value of `seed` sets the
#'                  seed for the random number generator.
#' @param num_processors An integer value. The value of `num_processors`
#'                  indicates the number of desired processors. Defaults to all
#'                  available cores using `parallelly::availableCores()`
#'
#' @returns A data frame object with three columns: the sequence name
#'   (`sequence_name`), the taxonomic assignment for each level (`taxonomy`),
#'   and the confidence associated with that classification level
#'   (`confidence`). The value in `confidence` represents the percentage of
#'   `num_bootstraps` that the classifier gave the same classification at that
#'   level.
#' @inherit build_kmer_database references
#' @export
#'
#' @examples
#' kmer_size <- 3
#' reference_names <- c("seq_1", "seq_2", "seq_3")
#' reference_sequences <- c("ATGCGCTA", "ATGCGCTC", "ATGCGCTC")
#' reference_genera <- c("A", "B", "B")
#' reference <- create_database(reference_names, reference_sequences, reference_genera)
#'
#' db <- build_kmer_database(reference, kmer_size)
#'
#' unknown_sequence <- c("ATGCGCTC", "ATGCGCTC", "ATGCGCTC")
#' unknown_id <- c("seq_1", "seq_2", "seq_3")
#'
#' classify_sequences(
#'   sequence = unknown_sequence,
#'   database = db,
#'   kmer_size = kmer_size,
#'   seed = 1,
#'   num_processors = 1
#' )
#'
#' classify_sequences(
#'   sequence = unknown_sequence,
#'   sequence_name = unknown_id,
#'   database = db,
#'   kmer_size = kmer_size
#' )
#'
#' classify_sequences(
#'   sequence = "ATGATCGAGTC",
#'   database = db,
#'   kmer_size = kmer_size,
#'   min_confidence = 95
#' )
#'
#' sequence_df <- data.frame(
#'   seq_names = unknown_id,
#'   seqs = unknown_sequence
#' )
#'
#' classify_sequences(sequence_df,
#'   sequence_name = seq_names, sequence = seqs,
#'   database = db, kmer_size = 3
#' )
#'

#' @export
#'@importFrom parallelly availableCores
classify_sequences <- function(
  data = NULL,
  sequence,
  sequence_name = as.character(seq_along(sequence)),
  database, kmer_size = 8, num_bootstraps = 100, min_confidence = 80,
  seed = 1,
  num_processors = parallelly::availableCores()
) {
  if (is.null(data)) {
    data <- FALSE
  }

  UseMethod("classify_sequences", data)
}



#' @export
classify_sequences.strollur <- function(
    data,
    sequence = NULL,
    sequence_name = NULL,
    database, kmer_size = 8, num_bootstraps = 100, min_confidence = 80,
    seed = 1,
    num_processors = parallelly::availableCores()) {
  if (!require_namespace("strollur")) {
    stop("Cannot use this function without strollur installed.")
  }
  if (strollur::xdev_has_sequence_taxonomy(data)) {
    warning(paste(
      "Your dataset already had taxonomy data, we are overwriting",
      "that with new classifications."
    ))
  }

  sequence <- strollur::xdev_get_sequences(data, degap = TRUE)
  sequence_name <- strollur::xdev_names(data, type = "sequence")
  df <- do_classify(
    sequence, sequence_name, database, kmer_size,
    num_bootstraps, min_confidence, seed, num_processors
  )

  if (!is.null(database$reference_data)) {
    strollur::xdev_add_references(
      data = data,
      table = database$reference_data,
      verbose = FALSE
    )
  }
  strollur::xdev_assign_sequence_taxonomy_tidy(
    data,
    table = df
  )
}


#' @export
classify_sequences.data.frame <- function(
    data, sequence = "sequence", sequence_name = "sequence_name",
    database, kmer_size = 8, num_bootstraps = 100, min_confidence = 80,
    seed = 1,
    num_processors = parallelly::availableCores()) {
  seq_colname <- as.character(substitute(sequence))
  sequence_name_colname <- as.character(substitute(sequence_name))
  df_name <- as.character(substitute(data))

  if (!(seq_colname %in% colnames(data))) {
    stop(paste0(
      "Column named '", seq_colname,
      "' does not exist in '", df_name, "'"
    ))
  }
  if (!(sequence_name_colname %in% colnames(data))) {
    stop(paste0(
      "Column named '", sequence_name_colname,
      "' does not exist in '", df_name, "'"
    ))
  }

  sequence <- data[, seq_colname]
  sequence_name <- data[, sequence_name_colname]
  do_classify(
    sequence, sequence_name, database, kmer_size,
    num_bootstraps, min_confidence, seed, num_processors
  )
}

#' @export
classify_sequences.default <- function(
  data = NULL,
  sequence,
  sequence_name = as.character(seq_along(sequence)),
  database, kmer_size = 8, num_bootstraps = 100, min_confidence = 80,
  seed = 1,
  num_processors = parallelly::availableCores()) {
  # are sequence and sequence_name the same length
  if (length(sequence) != length(sequence_name)) {
    message <- paste0(
      "You provided ", length(sequence), " sequence and ",
      length(sequence_name),
      " sequence_name `classify_sequences()` requires the ",
      "same number of sequence as sequence_name"
    )
    stop(message)
  }
  do_classify(
    sequence, sequence_name, database, kmer_size,
    num_bootstraps, min_confidence, seed, num_processors
  )
}

do_classify <- function(sequence, sequence_name, database, kmer_size,
                        num_bootstraps, min_confidence, seed,
                        num_processors) {
  df <- clean_sequence(sequence) |>
    # detect_kmers_rcpp(kmer_size) |>
    do_classify_rcpp(database$conditional_prob, num_bootstraps, kmer_size,
                     nrow(database$genera), seed, num_processors) |>
    consensus_bs_class(database$genera)


  df <- lapply(
    seq_along(along.with = df),
    \(x) {
      list(sequence_name = sequence_name[x], df[[x]])
    }
  ) |>
    lapply(X = _, as.data.frame) |>
    do.call(rbind, args = _)

  df <- df[df$confidence >= min_confidence, ]
  rownames(df) <- NULL

  df
}

require_namespace <- function(package_name) {
  requireNamespace("strollur", quietly = TRUE)
}
