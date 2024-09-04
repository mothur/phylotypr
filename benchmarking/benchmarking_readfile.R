file <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.fasta"

microbenchmark::microbenchmark(
  scan = scan(file, sep = "\n", what = character(), quiet = TRUE) |> stringi::stri_detect_regex("^>"),
  rl = readLines(file) |> stringi::stri_detect_regex("^>"),
  dt = data.table::fread(file, header = FALSE)$V1 |> stringi::stri_detect_regex("^>"),
  r_l = readr::read_lines(file) |> stringi::stri_detect_regex("^>"),
  v_l = vroom::vroom_lines(file) |> stringi::stri_detect_regex("^>"),
  times = 10
)
