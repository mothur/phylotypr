#include "detect_kmers.h"

inline uint64_t base_to_bits(char base) {
  switch (base) {
  case 'A': return 0b00;
  case 'C': return 0b01;
  case 'G': return 0b10;
  case 'T': return 0b11;
  default: return 0xFF; // Invalid character sentinel
  }
}

// Generates 2-bit hashes for all sliding k-mers in a sequence
std::vector<std::vector<uint64_t> > detect_kmers_rcpp(
    const std::vector<std::string>& sequence, const int kmer_size,
    long long start, long long end) {

  if(end == 0) { end = sequence.size(); }

  uint64_t max_kmer = std::pow(4, kmer_size) + 1;

  std::vector<std::vector<uint64_t> > kmers(end - start);
  int counter = 0;

  uint64_t mask = (kmer_size == 32) ? ~0ULL : (1ULL << (2 * kmer_size)) - 1;

  for(int s = start; s < end; s++){
    if (kmer_size > 32 || sequence[s].length() < kmer_size) return {};

    std::vector<uint64_t> hashes(sequence[s].length() - kmer_size + 1);

    uint64_t current_hash = 0;

    // Tracks the 0-based index of the most recent invalid character.
    // Initialized to -1 (represented via LL) meaning none seen yet.
    long long last_invalid_idx = -1;

    // Initialize the first k-mer window
    for (size_t i = 0; i < kmer_size; ++i) {
      uint64_t bits = base_to_bits(sequence[s][i]);
      if (bits == 0xFF) {
        last_invalid_idx = static_cast<long long>(i);
      }
      current_hash = (current_hash << 2) | (bits & 0b11); // Mask bits to safely shift
    }

    // Check if the first window contains any invalid character
    if (last_invalid_idx >= 0) {
      hashes[0] = max_kmer; //something we can remove later
    } else {
      hashes[0] = current_hash + 1;
    }

    int kmer_counter = 0;

    // Slide the window across the rest of the sequence
    for (size_t i = kmer_size; i < sequence[s].length(); ++i) {
      uint64_t bits = base_to_bits(sequence[s][i]);
      if (bits == 0xFF) {
        last_invalid_idx = static_cast<long long>(i);
      }

      current_hash = ((current_hash << 2) | (bits & 0b11)) & mask;

      // Calculate the starting boundary of the current k-mer window
      size_t window_start = i - kmer_size + 1;

      // If the last invalid character falls within the current window, it's invalid
      if (last_invalid_idx >= static_cast<long long>(window_start)) {
        hashes[++kmer_counter] = max_kmer; //something to remove later
      } else {
        hashes[++kmer_counter] = current_hash + 1;
      }
    }
    std::sort(hashes.begin(), hashes.end());
    auto hashes_iterator = std::unique(hashes.begin(), hashes.end());
    kmers[counter] = {hashes.begin(), hashes_iterator};

    if(kmers[counter][kmers[counter].size() - 1] == max_kmer){
      kmers[counter].pop_back();
    }
    counter++;
  }
  return kmers;
}
