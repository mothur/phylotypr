#include "classify_bs.h"

int classify_bs(const IntegerVector& unknown_kmers,
                const NumericMatrix& conditional_probs){

  int nrow = conditional_probs.nrow();
  int ncol = unknown_kmers.size();

  std::vector<double> probs(nrow, 0);

  for(int i = 0;i < ncol; ++i){

    const int col = unknown_kmers[i] - 1;

    for(int j = 0;j < nrow; ++j){
      probs[j] += *(conditional_probs.begin() + (nrow * col) + j);
    }

  }

  return(std::max_element(probs.begin(), probs.end())-probs.begin() + 1);
}
