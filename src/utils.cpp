#include "utils.h"

std::vector<pieceOfWork> Utils::divideWork(long long numItems, int& numProcessors) {

  // divide work between processors
  std::vector<pieceOfWork> work;

  if (numItems < numProcessors) { numProcessors = numItems; }
  size_t startIndex = 0;

  for (size_t remainingProcessors = numProcessors; remainingProcessors > 0;
       remainingProcessors--) {

    //case for last processor
    size_t numToProcess = numItems;
    if (remainingProcessors != 1) {
      numToProcess = ceil(numItems / remainingProcessors);
    }
    work.push_back(pieceOfWork(startIndex, (startIndex + numToProcess)));
    startIndex += numToProcess;
    numItems -= numToProcess;
  }
  return work;
}
