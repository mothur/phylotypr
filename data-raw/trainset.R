install_github("riffomonas/phylotypr")
library(phylotypr)

## code to prepare `trainset9_rdp` dataset goes here
temp_dir <- tempdir()

url <- "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.rdp.zip"
full_file_name <- basename(url)
temp_file_name <- paste0(temp_dir, "/", full_file_name)
download.file(url, temp_file_name)

unzip(temp_file_name, exdir = temp_dir)
fasta <- list.files(temp_dir,
                    pattern = ".rdp.fasta", full.names = TRUE)
taxonomy <- list.files(temp_dir,
                       pattern = ".rdp.tax", full.names = TRUE)

# fasta <- "benchmarking/trainset9_032012.rdp/trainset9_032012.rdp.fasta"
# taxonomy <- "benchmarking/trainset9_032012.rdp/trainset9_032012.rdp.tax"

fasta_df <- read_fasta(fasta)
genera <- read_taxonomy(taxonomy)

trainset9_rdp <- dplyr::inner_join(fasta_df, genera, by = "id")
trainset9_rdp <- trainset9_rdp[, c("id", "sequence", "taxonomy")]

usethis::use_data(trainset9_rdp, compress = "xz", overwrite = TRUE)




## code to prepare `trainset9_pds` dataset go here

temp_dir <- tempdir()

url <- "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip"
full_file_name <- basename(url)
temp_file_name <- paste0(temp_dir, "/", full_file_name)
download.file(url, temp_file_name)

unzip(temp_file_name, exdir = temp_dir)
fasta <- list.files(temp_dir,
                    pattern = ".pds.fasta", full.names = TRUE)
taxonomy <- list.files(temp_dir,
                       pattern = "pds.tax", full.names = TRUE)

# fasta <- "benchmarking/trainset9_032012.pds/trainset9_032012.pds.fasta"
# taxonomy <- "benchmarking/trainset9_032012.pds/trainset9_032012.pds.tax"

fasta_df <- read_fasta(fasta)
genera <- read_taxonomy(taxonomy)

trainset9_pds <- dplyr::inner_join(fasta_df, genera, by = "id")
trainset9_pds <- trainset9_pds[, c("id", "sequence", "taxonomy")]

usethis::use_data(trainset9_pds, compress = "xz", overwrite = TRUE)

