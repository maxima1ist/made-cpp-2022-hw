#pragma once

#include <fstream>
#include <utility>

namespace hw04 {

std::pair<size_t, size_t> get_count_good_and_bad(const std::string &filename,
                                                 const std::string &good_str,
                                                 const std::string &bad_str);

} // namespace hw04