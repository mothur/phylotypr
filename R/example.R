#' Get path to phylotypr example
#'
#' phylotypr comes bundled with some example files in its `inst/extdata`
#' directory. This function make them easy to access.
#'
#' @param path Name of file. If `NULL`, the example files will be listed.
#' @returns A string indicating path to the file listed in `path`. If `NULL` is
#'          given then the return value is a vector of file names in the
#'          `extdata/` directory
#' @export
#' @examples
#' phylotypr_example()
#' phylotypr_example("miseq_sop.fasta.gz")
phylotypr_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "phylotypr"))
  } else {
    system.file("extdata", path, package = "phylotypr", mustWork = TRUE)
  }
}
