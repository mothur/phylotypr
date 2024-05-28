#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calculate_log_probability(const NumericMatrix& kmer_genus_count,
                                        const NumericVector& word_specific_priors,
                                        const NumericVector& genus_counts){

  int n_kmers = kmer_genus_count.rows();
  int n_genera = kmer_genus_count.cols();

  NumericMatrix log_probability(n_kmers, n_genera);

  // log((kmer_genus_count + word_specific_priors) / (genus_counts + 1)

  for(int j = 0; j<n_genera; j++){
    int genus_counts_p1 = genus_counts[j] + 1;

    for(int i = 0; i<n_kmers; i++){
      log_probability(i, j) = log(
        (kmer_genus_count(i, j) + word_specific_priors(i))/(genus_counts_p1)
      );
    }
  }
  return log_probability;
}
