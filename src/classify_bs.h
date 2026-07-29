#ifndef CLASSIFY_BS_CPP
#define CLASSIFY_BS_CPP

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int classify_bs(const IntegerVector& unknown_kmers,
                const NumericMatrix& conditional_probs);

#endif
