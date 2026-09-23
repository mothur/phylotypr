#include <random> //test for marsenne twister
#include "rng.h"


std::vector<int> RNG::get_random_sample(const int low, const int high, const int amount) {
  std::vector<int> randomNumbers(amount);
  int range = high - low + 1;
  for (size_t i = 0; i < amount; i++) {
    randomNumbers[i] = low + (rng_engine() % range);
  }

  return randomNumbers;
}


