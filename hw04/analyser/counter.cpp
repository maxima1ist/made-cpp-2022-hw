#include "utils.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static constexpr const char *GOOD_STR = ":)";
static constexpr const char *BAD_STR = ":(";

std::string read_data_from_file(const std::string &filename) {
  std::ifstream fin(filename);
  std::string res((std::istreambuf_iterator<char>(fin)),
                  std::istreambuf_iterator<char>());
  fin.close();

  return res;
}

int main() {
  std::string filename;
  std::cin >> filename;

  const std::string chunk = read_data_from_file(filename);

  const auto good_count = hw04::count_substring_entries(
      {chunk.begin(), chunk.end()}, chunk.size(), GOOD_STR);
  const auto bad_count = hw04::count_substring_entries(
      {chunk.begin(), chunk.end()}, chunk.size(), BAD_STR);

  std::cout << good_count << ' ' << bad_count << '\n';

  return 0;
}