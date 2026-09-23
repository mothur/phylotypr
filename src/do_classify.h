#ifndef DO_CLASSIFY_CPP
#define DO_CLASSIFY_CPP

#include <Rcpp.h>
#include <vector>
#include <string>
#include "rng.h"

using namespace Rcpp;

struct classify_data{

  const std::vector<std::string>* sequence;
  long long start_index;
  long long end_index;
  const std::vector<double>* conditional_prob;
  const int num_bootstraps;
  const int kmer_size;
  const int n_genera;
  std::vector<std::vector<uint64_t> > bs_class;
  RNG random_number_generator;

  classify_data(const std::vector<std::string>& seqs,
                long long start_index,
                long long end_index,
                const std::vector<double>& con_prob,
                const int num_bootstraps,
                const int kmer_size,
                const int n_genera,
                const RNG rng):
     start_index(start_index), end_index(end_index),
     conditional_prob(&con_prob), num_bootstraps(num_bootstraps),
     kmer_size(kmer_size), n_genera(n_genera),
     random_number_generator(rng){
      bs_class.resize(end_index - start_index);
      sequence = &seqs;
    }

  ~classify_data() = default;
};


// [[Rcpp::export]]
std::vector<std::vector<uint64_t> > do_classify_rcpp(std::vector<std::string>& sequence,
                            const std::vector<double>& conditional_prob,
                            const int& num_bootstraps,
                            const int& kmer_size,
                            const int& n_genera,
                            const int& seed,
                            int num_processors = 1
                            );
#endif
