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

#include <exception>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "common/util/version.hpp"
#include "sql/ms_sql/ms_sql_util.hpp"
#include "sql/ms_sql/ms_sql_writer.hpp"

namespace {

void showUsage(const boost::program_options::options_description& desc) {
  std::cout << "Usage: top_sql_converter [options] spectrum-file-name"
            << std::endl;
  std::cout << desc << std::endl;
  std::cout << "Version: " << toppic::Version::getVersion() << std::endl;
}

}  // namespace

// Converts the MS1 peaks of an mzML/mzXML file into an SQLite database for
// the 3D visualization (see sql/ms_sql/ms_sql_writer.hpp for the layout).
int main(int argc, char* argv[]) {
  namespace po = boost::program_options;
  double mz_size = toppic::MsSqlWriter::DEFAULT_MZ_SIZE;
  double rt_divider = toppic::MsSqlWriter::DEFAULT_RT_DIVIDER;
  std::string spec_file_name;

  po::options_description visible_desc("Options");
  visible_desc.add_options()("help,h", "Print this help message.")(
      "mz-size,m", po::value<double>(&mz_size),
      "<a positive number>. Set the m/z width of a grid block in the first "
      "down-sampled peak layer. The default value is 0.05 m/z.")(
      "rt-divider,r", po::value<double>(&rt_divider),
      "<a positive number>. The retention time height of a grid block in the "
      "first down-sampled peak layer is the average MS1 scan interval divided "
      "by this value. The default value is 1.");
  po::options_description hidden_desc("Hidden options");
  hidden_desc.add_options()("spectrum-file-name",
                            po::value<std::string>(&spec_file_name)->required(),
                            "Spectrum file name.");
  po::options_description desc;
  desc.add(visible_desc).add(hidden_desc);
  po::positional_options_description positional;
  positional.add("spectrum-file-name", 1);

  try {
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
                  .options(desc)
                  .positional(positional)
                  .run(),
              vm);
    if (vm.count("help")) {
      showUsage(visible_desc);
      return 0;
    }
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    showUsage(visible_desc);
    return 1;
  }
  if (mz_size <= 0.0 || rt_divider <= 0.0) {
    std::cerr << "ERROR: mz-size and rt-divider must be positive." << std::endl;
    return 1;
  }

  std::cout << "top_sql_converter " << toppic::Version::getVersion()
            << std::endl;
  try {
    toppic::ms_sql_util::convert(spec_file_name, mz_size, rt_divider);
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
