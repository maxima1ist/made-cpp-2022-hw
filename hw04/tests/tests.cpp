#include "analyser.hpp"
#include "utils.hpp"

#include <gtest/gtest.h>

#include <iostream>

static constexpr const char *GOOD_STR = ":)";
static constexpr const char *BAD_STR = ":(";

namespace {

std::vector<char> get_vec_from_str(const std::string &str) {
  return {str.begin(), str.end()};
}

} // namespace

TEST(unit_test, count_substring_entries) {
  std::string chunk;
  ASSERT_EQ(
      hw04::count_substring_entries(get_vec_from_str(chunk), chunk.size(), ""),
      0);

  chunk = "Some text";
  ASSERT_EQ(
      hw04::count_substring_entries(get_vec_from_str(chunk), chunk.size(), ""),
      0);

  ASSERT_EQ(hw04::count_substring_entries(get_vec_from_str(chunk), chunk.size(),
                                          "text"),
            1);

  chunk = "Some long long long long long long long long long long long long "
          "long long long long long long long long long long long long long "
          "long long long long long long long long long long long long long "
          "long long long long long long long long long long long long long "
          "long long long long long long long long long long long long long "
          "long long long long long long long long line";
  ASSERT_EQ(hw04::count_substring_entries(get_vec_from_str(chunk), chunk.size(),
                                          "long"),
            72);
}

TEST(unit_test, get_count_good_and_bad) {
  auto lhs = hw04::get_count_good_and_bad("", GOOD_STR, BAD_STR);
  ASSERT_EQ(lhs.first, 0);
  ASSERT_EQ(lhs.second, 0);

  lhs = hw04::get_count_good_and_bad("../../data/sample_1mb.txt", GOOD_STR,
                                     BAD_STR);
  ASSERT_EQ(lhs.first, 8320);
  ASSERT_EQ(lhs.second, 4160);

  lhs = hw04::get_count_good_and_bad("../../data/sample_10mb.txt", GOOD_STR,
                                     BAD_STR);
  ASSERT_EQ(lhs.first, 27115);
  ASSERT_EQ(lhs.second, 15840);

  lhs = hw04::get_count_good_and_bad("../../data/sample_100mb.txt", GOOD_STR,
                                     BAD_STR);
  ASSERT_EQ(lhs.first, 325380);
  ASSERT_EQ(lhs.second, 190080);
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}