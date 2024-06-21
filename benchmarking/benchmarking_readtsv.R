taxonomy <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.tax"

genera <- readr::read_tsv(taxonomy,
                          col_names = c("accession", "taxonomy")) |>
  dplyr::mutate(taxonomy = stringi::stri_replace_all_regex(taxonomy, ";$", ""))


microbenchmark::microbenchmark(
  r.d = read.delim(taxonomy, header = FALSE,
             col.names = c("accession", "taxonomy"),
             colClasses = c("character", "character")) |>
    dplyr::mutate(taxonomy = stringi::stri_replace_last_regex(taxonomy, ";$", "")),
  r_t = readr::read_tsv(taxonomy,
                  col_names = c("accession", "taxonomy"),
                  col_types = cols(.default = col_character())) |>
    dplyr::mutate(taxonomy = stringi::stri_replace_last_regex(taxonomy, ";$", "")),
  v = vroom::vroom(taxonomy, delim = "\t",
               col_names = c("accession", "taxonomy"),
               col_types = cols(.default = col_character())) |>
    dplyr::mutate(taxonomy = stringi::stri_replace_last_regex(taxonomy, ";$", "")),
  d.t = data.table::fread(taxonomy, sep = "\t", header = FALSE,
                    col.names = c("accession", "taxonomy")) |>
    dplyr::mutate(taxonomy = stringi::stri_replace_last_regex(taxonomy, ";$", ""))

)
