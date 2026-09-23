#ifndef UTILS_CPP
#define UTILS_CPP

#include <Rcpp.h>
#include <vector>
#include <string>

struct pieceOfWork {
  long long start;
  long long end;
  pieceOfWork(long long i, long long j) : start(i), end(j) {}
  pieceOfWork() { start = 0; end = 0; }
  ~pieceOfWork() {}
};

class Utils {
public:
  static std::vector<pieceOfWork> divideWork(long long numItems, int& numProcessors);

};

#endif
