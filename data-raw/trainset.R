# install_github("mothur/phylotypr")
library(phylotypr)



## code to prepare `rdp_trainset9` dataset goes here
temp_dir <- tempdir()

url <- "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.rdp.zip"
full_file_name <- basename(url)
temp_file_name <- paste0(temp_dir, "/", full_file_name)
download.file(url, temp_file_name)

unzip(temp_file_name, exdir = temp_dir)
fasta <- list.files(temp_dir,
  pattern = ".rdp.fasta", full.names = TRUE
)
taxonomy <- list.files(temp_dir,
  pattern = ".rdp.tax", full.names = TRUE
)

fasta_df <- read_fasta(fasta)
genera <- read_taxonomy(taxonomy)

rdp_trainset9 <- dplyr::inner_join(fasta_df, genera, by = "sequence_name")
rdp_trainset9 <- rdp_trainset9[, c("sequence_name", "sequence", "taxonomy")]

reference_data <-
  data.frame(
    name = "rdp_trainset9", # rdp_trainset_10
    vendor = "rdp", # rdp, pds, silva, greengenes, ncbi
    version = "9", # 10_pds, v8,
    usage = "classification", # "classification"
    note = paste0("phylotypr v.", packageVersion("phylotypr")), # phylotyprrefdata version
    method_url = "https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/", # "mothur blog post creating reference?"
    documentation_url = "https://mothur.org/wiki/rdp_reference_files/#version-9, ",
    citation =
      "@article{doi:10.1128/AEM.00062-07,
                    author = {Qiong Wang and George M. Garrity and James M. Tiedje and James R. Cole},
                    title = {Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy},
                    journal = {Applied and Environmental Microbiology},
                    volume = {73},
                    number = {16},
                    pages = {5261-5267},
                    year = {2007},
                    doi = {10.1128/AEM.00062-07},

                    URL = {https://journals.asm.org/doi/abs/10.1128/aem.00062-07},
                    eprint = {https://journals.asm.org/doi/pdf/10.1128/aem.00062-07}
                  }"
  )

rdp_trainset9 <- create_database(sequence_name = rdp_trainset9$sequence_name,
                                 sequence = rdp_trainset9$sequence,
                                 taxonomy = rdp_trainset9$taxonomy,
                                 reference_data = reference_data)

usethis::use_data(rdp_trainset9, compress = "xz", overwrite = TRUE)


## code to prepare `trainset9_pds` dataset go here

temp_dir <- tempdir()

url <- "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip"
full_file_name <- basename(url)
temp_file_name <- paste0(temp_dir, "/", full_file_name)
download.file(url, temp_file_name)

unzip(temp_file_name, exdir = temp_dir)
fasta <- list.files(temp_dir,
  pattern = ".pds.fasta", full.names = TRUE
)
taxonomy <- list.files(temp_dir,
  pattern = "pds.tax", full.names = TRUE
)

fasta_df <- read_fasta(fasta)
genera <- read_taxonomy(taxonomy)

rdp_trainset9_pds <- dplyr::inner_join(fasta_df, genera, by = "sequence_name")
rdp_trainset9_pds <- rdp_trainset9_pds[, c("sequence_name", "sequence", "taxonomy")]


reference_data <-
  data.frame(
    name = "rdp_trainset9_pds", # rdp_trainset_10
    vendor = "rdp", # rdp, pds, silva, greengenes, ncbi
    version = "9_pds", # 10_pds, v8,
    usage = "classification", # "classification"
    note = paste0("phylotypr v.", packageVersion("phylotypr")), # phylotyprrefdata version
    method_url = "https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/", # "mothur blog post creating reference?"
    documentation_url = "https://mothur.org/wiki/rdp_reference_files/#version-9, ",
    citation =
      "@article{doi:10.1128/AEM.00062-07,
                    author = {Qiong Wang and George M. Garrity and James M. Tiedje and James R. Cole},
                    title = {Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy},
                    journal = {Applied and Environmental Microbiology},
                    volume = {73},
                    number = {16},
                    pages = {5261-5267},
                    year = {2007},
                    doi = {10.1128/AEM.00062-07},

                    URL = {https://journals.asm.org/doi/abs/10.1128/aem.00062-07},
                    eprint = {https://journals.asm.org/doi/pdf/10.1128/aem.00062-07}
                  }"
  )


rdp_trainset9_pds <- create_database(
  sequence_name = rdp_trainset9_pds$sequence_name,
  sequence = rdp_trainset9_pds$sequence,
  taxonomy = rdp_trainset9_pds$taxonomy,
  reference_data = reference_data)

usethis::use_data(rdp_trainset9_pds, compress = "xz", overwrite = TRUE)
