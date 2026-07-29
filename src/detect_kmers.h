#ifndef DETECT_KMERS_CPP
#define DETECT_KMERS_CPP

#include <Rcpp.h>
#include <string>
#include <vector>
#include <iostream>
#include <string_view>
#include <optional>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
std::vector<std::vector<uint64_t> > DetectKmers(std::vector<std::string>& sequence,
                                                const int kmer_size = 8);


#endif
