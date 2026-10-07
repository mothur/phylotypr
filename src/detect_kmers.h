#ifndef DETECT_KMERS_CPP
#define DETECT_KMERS_CPP

#include <Rcpp.h>
#include <string>
#include <vector>
#include <iostream>
#include <string_view>
#include <optional>

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::vector<uint64_t> > detect_kmers_rcpp(
    const std::vector<std::string>& sequence,
    const int kmer_size = 8,
    long long start = 0, long long end = 0);


#endif
