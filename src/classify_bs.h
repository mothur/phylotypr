#ifndef CLASSIFY_BS_CPP
#define CLASSIFY_BS_CPP

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int classify_bs_rcpp(const std::vector<int>& unknown_kmers,
                     const std::vector<double>& conditional_probs,
                     const int& n_genera);

#endif
