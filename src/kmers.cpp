#include <Rcpp.h>
#include <chrono>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> calculate_log_probability_rcpp(// produces probability of seeing kmer in genus
    const Rcpp::List& detect_list, // index of kmers detected in each sequence
    const NumericVector& word_specific_priors, // probability of seeing kmer across all sequences
    const NumericVector& genera, // the genus index for each sequence
    const NumericVector& genus_counts){ //number of sequences in each genus

  int n_genera = genus_counts.size();
  int n_sequences = detect_list.size();
  int n_kmers = word_specific_priors.size();

  std::vector<std::vector<int> > kmer_genus_count(n_kmers, std::vector<int>(n_genera, 0));

  for(int i = 0; i < n_sequences; i++){
    const int genera_index = genera[i] - 1;
    const std::vector<int>& sequence_kmers = Rcpp::as<std::vector<int>>(detect_list[i]);
    for(int j = 0; j <  sequence_kmers.size(); j++){
      const int sequence_kmer_index =  sequence_kmers[j] - 1;
      kmer_genus_count[sequence_kmer_index][genera_index]++; // need to -1 because 1 vs 0 index; probably a problem elsewhere
    }
  }

  std::vector<double> log_probability(n_kmers * n_genera);

  for(int i = 0; i < n_kmers; i++){
    int i_n_genera = i * n_genera;
    const std::vector<int>& current_kmer_count = kmer_genus_count[i];
    const double word_prior = word_specific_priors(i);
    for(int j = 0; j<n_genera; j++){
      const int genus_counts_p1 = genus_counts[j] + 1;
      log_probability[i_n_genera + j] = std::log(
        (current_kmer_count[j] + word_prior)/(genus_counts_p1)
      );

    }
  }

  return log_probability;
}
