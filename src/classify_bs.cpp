#include "classify_bs.h"

int classify_bs_rcpp(const std::vector<int>& unknown_kmers,
                     const std::vector<double>& conditional_probs,
                     const int& n_genera){

  int nrow = n_genera;
  int ncol = unknown_kmers.size();

  std::vector<double> probs(nrow, 0);
  for(int i = 0;i < ncol; ++i){
    const int nrow_col = nrow * (unknown_kmers[i] - 1);
    for(int j = 0;j < nrow; ++j){
      probs[j] += conditional_probs[nrow_col + j];
    }
  }
  return(std::max_element(probs.begin(), probs.end())-probs.begin() + 1);
}
