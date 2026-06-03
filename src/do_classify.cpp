#ifndef do_classify_cpp
#define do_classify_cpp

#include <Rcpp.h>
#include <vector>
#include <string>
#include "classify_bs.h"

// [[Rcpp::export]]
Rcpp::List do_classify_rcpp(
    const std::vector<std::vector<int> >& sequence_kmers,
    const Rcpp::NumericMatrix& conditional_prob,
    const int& num_bootstraps,
    const int& kmer_size) {

  Rcpp::List bs_class(sequence_kmers.size());
  int sequence_index = 0;

  for(const std::vector<int>& sequence:sequence_kmers){

    int sample_n_kmers = floor(sequence.size() / kmer_size);
    Rcpp::IntegerVector rcpp_sequence = Rcpp::wrap(sequence);
    Rcpp::IntegerVector result(num_bootstraps);

    for(int i = 0;i<num_bootstraps; i++){
      Rcpp::IntegerVector bs_kmers = Rcpp::sample(rcpp_sequence, sample_n_kmers,
                                                  true);

      result[i] = classify_bs(bs_kmers, conditional_prob);
    }
    bs_class[sequence_index] = result;
    sequence_index++;
  }

  return bs_class;
}

#endif
