// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// topzstd: compresses the SQLite databases of TopFD, TopDIA, toppic and
// topmg (sqlite/zstd_db; toppic_test experiment 172), or with -d
// decompresses them. Those tools compress their database by default; the
// tool converts the databases of earlier versions, or of runs with
// --disable-zstd, and gives the rows back for programs that read SQLite
// without the sqlite-zstd extension.

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "common/util/logger.hpp"
#include "common/util/version.hpp"
#include "sql/zstd_db.hpp"

namespace {

void showUsage(const boost::program_options::options_description& desc) {
  std::cout << "Usage: topzstd [options] database-file-name ..." << std::endl;
  std::cout << desc << std::endl;
  std::cout << "Version: " << toppic::Version::getVersion() << std::endl;
}

}  // namespace

int main(int argc, char* argv[]) {
  toppic::logger::log_level = 5;
  namespace po = boost::program_options;
  bool decompress = false;
  int level = 3;
  std::string output;
  std::vector<std::string> inputs;
  po::options_description display_desc("Options");
  display_desc.add_options()("help,h", "Print the help message.")(
      "decompress,d",
      "Decompress: write the database with its peak and envelope tables back "
      "as rows, as it was before compression.")(
      "level,l", po::value<int>(&level),
      "<1-22>. The zstd compression level. Default value: 3.")(
      "output,o", po::value<std::string>(&output),
      "<a file name>. Write the result to this file instead of replacing "
      "the input file (one input file only).");
  po::options_description desc("Options");
  desc.add(display_desc);
  desc.add_options()("database-file-name",
                     po::value<std::vector<std::string>>(&inputs));
  po::positional_options_description positional;
  positional.add("database-file-name", -1);
  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv)
                  .options(desc)
                  .positional(positional)
                  .run(),
              vm);
    if (vm.count("help")) {
      showUsage(display_desc);
      return 0;
    }
    po::notify(vm);
  } catch (po::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    showUsage(display_desc);
    return 1;
  }
  decompress = vm.count("decompress") > 0;
  if (inputs.empty()) {
    std::cerr << "ERROR: no database file name." << std::endl << std::endl;
    showUsage(display_desc);
    return 1;
  }
  if (level < 1 || level > 22) {
    std::cerr << "ERROR: the zstd level must be 1 - 22." << std::endl;
    return 1;
  }
  if (!output.empty() && inputs.size() > 1) {
    std::cerr << "ERROR: --output takes one input file." << std::endl;
    return 1;
  }
  // every input is checked before any is changed
  for (const std::string& name : inputs) {
    if (!std::filesystem::exists(name)) {
      std::cerr << "ERROR: " << name << " does not exist." << std::endl;
      return 1;
    }
    if (!toppic::zstd_db::isSqliteFile(name)) {
      std::cerr << "ERROR: " << name << " is not an SQLite database."
                << std::endl;
      return 1;
    }
    bool compressed = toppic::zstd_db::isCompressed(name);
    if (!decompress && compressed) {
      std::cerr << "ERROR: " << name << " is compressed already." << std::endl;
      return 1;
    }
    if (decompress && !compressed) {
      std::cerr << "ERROR: " << name << " is not compressed." << std::endl;
      return 1;
    }
  }
  for (const std::string& name : inputs) {
    std::string dest = output.empty() ? name : output;
    std::cout << (decompress ? "Decompressing " : "Compressing ") << name
              << (dest == name ? "" : " to " + dest) << " started."
              << std::endl;
    toppic::zstd_db::Stat stat =
        decompress ? toppic::zstd_db::decompressDb(name, dest)
                   : toppic::zstd_db::compressDb(name, dest, level);
    std::cout << (decompress ? "Decompressing " : "Compressing ") << name
              << " finished: " << std::fixed << std::setprecision(1)
              << stat.bytes_before / 1e6 << " MB -> " << stat.bytes_after / 1e6
              << " MB (" << std::setprecision(1)
              << 100.0 * stat.bytes_after /
                     (stat.bytes_before > 0 ? stat.bytes_before : 1)
              << "%), " << stat.seconds << " s." << std::defaultfloat
              << std::endl;
  }
  return 0;
}
