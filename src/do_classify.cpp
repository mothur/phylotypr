#include "do_classify.h"
#include "detect_kmers.h"
#include "classify_bs.h"
#include "utils.h"
#include <RcppThread.h>

void driver_classify(classify_data* data){

  const std::vector<std::vector<uint64_t> > sequence_kmers = detect_kmers_rcpp(
    *data->sequence, data->kmer_size, data->start_index, data->end_index);

  std::vector<std::vector<int> > bs_class(sequence_kmers.size());
  int sequence_index = 0;

  for(const std::vector<uint64_t>& sequence:sequence_kmers){

    const int sequence_size = sequence.size();
    int sample_n_kmers = floor(sequence_size / data->kmer_size);
    // std::vector<uint64_t> result(data->num_bootstraps);
    data->bs_class[sequence_index].resize(data->num_bootstraps);

    for(int i = 0;i<data->num_bootstraps; i++){
      std::vector<int> random_indices = data->random_number_generator.get_random_sample(0,  sequence_size - 1, sample_n_kmers);
      std::vector<int> bs_kmers(random_indices.size());

      for(int j = 0; j< bs_kmers.size(); j++){
        bs_kmers[j] = static_cast<int>(sequence[random_indices[j]]);
      }

      //should return indices into the kmers
      data->bs_class[sequence_index][i] = classify_bs_rcpp(bs_kmers,
                                   *data->conditional_prob, data->n_genera);
    }

    // data->bs_class[sequence_index] = result;
    sequence_index++;
  }

}

//[Depends(RcppThread)]
std::vector<std::vector<uint64_t> > do_classify_rcpp(
    std::vector<std::string>& sequence,
    const std::vector<double>& conditional_prob,
    const int& num_bootstraps,
    const int& kmer_size,
    const int& n_genera,
    const int& seed,
    int num_processors)   {

  std::vector<pieceOfWork> work = Utils::divideWork(sequence.size(),
                                                    num_processors);
  std::vector<RcppThread::Thread*> workers;
  std::vector<classify_data*> data;

  for(int i = 1; i < num_processors; i++){
    classify_data* data_bundle = new classify_data(
      sequence, work[i].start, work[i].end, conditional_prob,
      num_bootstraps, kmer_size, n_genera, RNG(seed + i)
    );

    data.push_back(data_bundle);
    workers.push_back(new RcppThread::Thread(driver_classify, data_bundle));
  }

  classify_data* data_bundle = new classify_data(
    sequence, work[0].start, work[0].end, conditional_prob,
    num_bootstraps, kmer_size, n_genera, RNG(seed)
  );

  driver_classify(data_bundle);

// need to pull back values from data_bundles and data
  std::vector<std::vector<uint64_t> > results = data_bundle->bs_class;

  for(int i = 0; i < num_processors - 1; i++){

    workers[i]->join();

    results.insert(results.end(),
                   data[i]->bs_class.begin(), data[i]->bs_class.end());

    delete workers[i];
    delete data[i];
  }

  delete data_bundle;

  return results;
}
