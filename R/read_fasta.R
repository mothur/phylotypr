#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
read_fasta <- function(file) {

  fasta_data <- scan(file,
                     sep = "\n",
                     what = character(), quiet = TRUE)

  n_lines <- length(fasta_data)

  id <- character(n_lines)
  sequence <- character(n_lines)
  comment <- character(n_lines)

  index <- 1
  temp_sequence <- ""

  for(line in seq_along(fasta_data)){

    if(stringi::stri_detect_regex(fasta_data[line], "^>")){
      id_line <- fasta_data[line]

      id[index] <- id_line |>
        stringi::stri_replace_all_regex("^>(\\w*)\\s*.*", "$1")

      comment[index] <- id_line |>
        stringi::stri_replace_all_regex("^>\\w*\\s*(.*)", "$1")

      sequence[index] <- paste(temp_sequence, collapse = "")

      temp_sequence <- ""
      index <- index + 1

    } else {
      temp_sequence <- c(temp_sequence, fasta_data[line])
    }


  }

  sequence[index] <- paste(temp_sequence, collapse = "")

  data.frame(id = id[1:(index - 1)],
             sequence = sequence[2: index],
             comment = comment[1:(index - 1)])
}
