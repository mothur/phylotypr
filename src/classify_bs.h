#ifndef classify_bs_cpp
#define classify_bs_cpp

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int classify_bs(const IntegerVector& unknown_kmers,
                const NumericMatrix& conditional_probs);

#endif
