#include "analyser.hpp"

#include "utils.hpp"

#include <ios>

namespace hw04 {

std::pair<size_t, size_t> get_count_good_and_bad(const std::string &filename,
                                                 const std::string &good_str,
                                                 const std::string &bad_str) {
  std::ifstream fin(filename);

  size_t good_count = 0;
  size_t bad_count = 0;
  if (!fin.is_open()) {
    return std::make_pair(good_count, bad_count);
  }

  std::vector<char> chunk(4 * 1024 * 1024, 0);
  char last_char;
  while (!fin.eof()) {
    fin.read(chunk.data(), chunk.size());
    std::streamsize chunk_size = fin.gcount();

    if (last_char == ':') {
      if (chunk[0] == ')') {
        ++good_count;
      } else if (chunk[0] == '(') {
        ++bad_count;
      }
    }
    last_char = chunk[chunk_size - 1];

    good_count += count_substring_entries(chunk, chunk_size, good_str);
    bad_count += count_substring_entries(chunk, chunk_size, bad_str);
  }
  fin.close();

  return std::make_pair(good_count, bad_count);
}

} // namespace hw04