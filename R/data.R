#' RDP training set v9
#'
#' The sequence and taxonomy data for the 10,049 sequences found in the
#' Ribosomal Database Project's trainset9_032012 training set for use with the
#' naive Bayesian classifier as implemented in the `{phylyotypr}` R package.
#' Originally released by the RDP in September 2012. The `rdp` version contains
#' the same sequences as provided by the official RDP version (9,665 bacterial
#' and 384 archaeal). The `pds` version contains extra eukaryotic sequences
#' including 119 chloroplasts and mitochondria (10,168 total sequences). See the
#' mothur reference file page in "Sources" for more information. Be sure to see
#' the mothur GitHub project where you can find the phylotyprrefdata package
#' (https://github.com/mothur/phylotyprrefdata) for access to other taxonomic
#' reference data.
#'
#'
#' @format
#' A `classification_database` object that holds two slots: seq_taxonomy_data
#' and `reference_data.` The seq_taxonomy_data slot is a data.frame containing
#' \describe{
#' \item{sequence_name}{Sequence accession identifier}
#' \item{sequence}{DNA sequence string}
#' \item{taxonomy}{Taxonomic string with each level separated with a `;`}
#' }
#' The `reference_data` slot contains metadata about the database and includes
#' fields:
#' \describe{
#' \item{reference_vendor}{Entity that created original database (e.g., RDP)}
#' \item{reference_versions}{The version of database provided by original
#'  entity}
#' \item{reference_usages}{For classification or alignment}
#' \item{reference_notes}{Any additional notes about the database. Currently the
#'  version of phylotypr that provided the database}
#' \item{reference_urls}{URLs where database can be obtained}
#' \item{reference_citation}{Citation for database in BibTeX format}
#' }
#' @source
#' * [mothur-formatted files](https://mothur.org/wiki/rdp_reference_files/)
#' * [RDP sourceforge page](https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo9_rawtrainingdata.zip/download) # nolint: line_length_linter
#'
"rdp_trainset9"


#' @rdname rdp_trainset9
"rdp_trainset9_pds"
