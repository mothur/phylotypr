#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame df_rcpp(int n) {

  Vector kmer = sample(64, n, true);
  Vector genus = sample(15, n, true);
  Vector count = sample(10, n, true);

  DataFrame df = DataFrame::create(Named("kmer") = kmer,
                                   Named("genus") = genus,
                                   Named("count") = count);

  return df;
}
