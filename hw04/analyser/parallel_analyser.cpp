#include "analyser.hpp"

#include "config.h"
#include "utils.hpp"

#include <boost/process.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <unordered_map>

namespace fs = std::filesystem;
namespace bp = boost::process;

namespace {

struct ChildProcess {
  bp::child ps;
  bp::ipstream out;
  std::string filename;
};

void write_data_to_file(const std::string &filename, const char *data,
                        size_t data_size) {
  std::ofstream fout(filename);
  fout.write(data, data_size);
  fout.close();
}

} // namespace

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

  fs::path file_path(filename);
  auto file_size = std::filesystem::file_size(file_path);

  std::unordered_map<size_t, ChildProcess> child_processes;
  size_t child_id = 0;
  boost::asio::io_context io;
  auto work = make_work_guard(io);

  const size_t ps_count = std::thread::hardware_concurrency();
  std::vector<char> chunk(file_size / ps_count, 0);
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

    child_processes[child_id].filename =
        ANALYZER_PATH "/filetmp" + std::to_string(child_id) + ".txt";
    write_data_to_file(child_processes[child_id].filename,
                       reinterpret_cast<char *>(chunk.data()), chunk_size);

    const std::string exec_cmd =
        ANALYZER_PATH "/counter " + child_processes[child_id].filename;

    bp::pstream p;
    p << child_processes[child_id].filename << std::endl;
    auto &out = child_processes[child_id].out;
    auto &child_ps = child_processes[child_id].ps;
    child_ps = bp::child(exec_cmd,
                         bp::on_exit([&out, &good_count, &bad_count](
                                         int, const std::error_code &) {
                           std::string good_count_str, bad_count_str;
                           out >> good_count_str >> bad_count_str;

                           good_count += std::stoi(good_count_str);
                           bad_count += std::stoi(bad_count_str);
                         }),
                         io, bp::std_out > out, bp::std_in = p);

    ++child_id;
  }
  fin.close();

  work.reset();

  while (io.run_one()) {
  }

  for (auto &child_process : child_processes) {
    std::remove(child_process.second.filename.c_str());
  }

  return std::make_pair(good_count, bad_count);
}

} // namespace hw04