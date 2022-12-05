#include "utils.hpp"

namespace hw04 {

size_t count_substring_entries(const std::vector<char> &chunk,
                               size_t chunk_size, const std::string &substr) {
  size_t count = 0;

  if (chunk.empty() || substr.empty()) {
    return count;
  }

  for (size_t i = 0; i < chunk_size; ++i) {
    size_t pos = 0;
    while (pos < substr.size() && chunk[i] == substr[pos]) {
      ++pos;
      if (pos != substr.size()) {
        ++i;
      }
    }

    if (pos == substr.size()) {
      ++count;
    }
  }

  return count;
}

}