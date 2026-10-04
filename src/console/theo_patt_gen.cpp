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

// theo_patt_gen writes TopFD's averagine table res/base_data/theo_patt.txt
// (common/emass/averagine_table). A developer tool: not built by default
// (make theo_patt_gen) and not installed.

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "common/emass/averagine_table.hpp"
#include "common/emass/isotope_data.hpp"
#include "common/util/version.hpp"

namespace {

// The largest value of --sulphur-num.
const int kMaxSulphurNum = 5;

}  // namespace

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;
  std::string output = "theo_patt.txt";
  std::string composition = "human";
  std::string sulphur_str = "0";
  po::options_description desc("Options");
  desc.add_options()("help,h", "Print the help message.")(
      "output,o", po::value<std::string>(&output),
      "<file name>. The averagine table. Default: theo_patt.txt.")(
      "composition,c", po::value<std::string>(&composition),
      "<human|human_no_sulphur|senko>. The averagine: human, the mean "
      "residue of the UniProt human reference proteome (UP000005640, "
      "release 2026_03); human_no_sulphur, the same without cysteine and "
      "methionine (C, H, N and O only); senko, Senko, Beu and McLafferty's "
      "(JASMS 6:229, 1995). Default: human.")(
      "sulphur-num", po::value<std::string>(&sulphur_str),
      "<0-5>. The number k of sulphur atoms in every formula: the "
      "composition is scaled to the mass minus k 32S masses and k S are "
      "added. Only with human_no_sulphur when k > 0. Default: 0.");
  po::variables_map vm;
  try {
    // no positional arguments: a stray one is an error, not ignored
    po::positional_options_description no_positional;
    po::store(po::command_line_parser(argc, argv)
                  .options(desc)
                  .positional(no_positional)
                  .run(),
              vm);
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (vm.count("help")) {
    std::cout << "Usage: theo_patt_gen [options]" << std::endl;
    std::cout << "Writes TopFD's averagine table (res/base_data/theo_patt.txt):"
              << std::endl
              << "11,000 entries, one per 10 Da, by the emass method with "
                 "AME2020 isotope"
              << std::endl
              << "masses and IUPAC representative abundances." << std::endl;
    std::cout << desc << std::endl;
    std::cout << "Version: " << toppic::Version::getVersion() << std::endl;
    return EXIT_SUCCESS;
  }
  toppic::emass::Composition comp;
  if (composition == "human") {
    comp = toppic::emass::kHumanAveragine;
  } else if (composition == "human_no_sulphur") {
    comp = toppic::emass::kHumanNoSulphurAveragine;
  } else if (composition == "senko") {
    comp = toppic::emass::kSenkoAveragine;
  } else {
    std::cerr << "Error: unknown composition " << composition
              << "; use human, human_no_sulphur or senko." << std::endl;
    return EXIT_FAILURE;
  }
  if (sulphur_str.empty() ||
      !std::all_of(sulphur_str.begin(), sulphur_str.end(),
                   [](unsigned char c) { return std::isdigit(c); })) {
    std::cerr << "Error: --sulphur-num " << sulphur_str
              << " is not a non-negative integer." << std::endl;
    return EXIT_FAILURE;
  }
  if (sulphur_str.size() > 2 || std::stoi(sulphur_str) > kMaxSulphurNum) {
    std::cerr << "Error: --sulphur-num " << sulphur_str << " is above "
              << kMaxSulphurNum << ", the largest value." << std::endl;
    return EXIT_FAILURE;
  }
  int sulphur_num = std::stoi(sulphur_str);
  if (sulphur_num > 0 && composition != "human_no_sulphur") {
    std::cerr << "Error: --sulphur-num " << sulphur_num
              << " needs the human_no_sulphur composition (-c "
                 "human_no_sulphur)." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Writing the averagine table (" << composition
            << " averagine";
  if (sulphur_num > 0) {
    std::cout << ", " << sulphur_num << " S per formula";
  }
  std::cout << ") to " << output << " - started." << std::endl;
  toppic::emass::writeAveragineTable(output, comp, sulphur_num);
  std::cout << "Writing the averagine table - finished." << std::endl;
  return EXIT_SUCCESS;
}
