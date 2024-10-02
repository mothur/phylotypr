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
#' A data frame with 3 columns. Each row represents a different
#' sequence:
#' \describe{
#'   \item{id}{Sequence accession identifier}
#'   \item{sequence}{DNA sequence string}
#'   \item{taxonomy}{Taxonomic string with each level separated with a `;`}
#' }
#'
#' @source
#' * [mothur-formatted files](https://mothur.org/wiki/rdp_reference_files/)
#' * [RDP sourceforge page](https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo9_rawtrainingdata.zip/download) # nolint: line_length_linter
#'
"trainset9_rdp"


#' @rdname trainset9_rdp
"trainset9_pds"
