#ifndef RNG_H
#define RNG_H
#include <vector>
#include <random>

class RNG {
public:
  RNG(int seed):rng_engine(std::mt19937(seed)){}
  std::vector<int> get_random_sample(const int low, const int high,
                                  const int amount);
private:
  std::mt19937 rng_engine;

};


#endif //RNG_H
