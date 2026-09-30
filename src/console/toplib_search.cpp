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

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "common/util/logger.hpp"
#include "console/toplib_search_argument.hpp"
#include "toplib/lib_msalign_reader.hpp"
#include "toplib/lib_numeric.hpp"
#include "toplib/lib_search.hpp"
#include "toplib/lib_spectrum.hpp"

using namespace toppic;

int main(int argc, char* argv[]) {
  toppic::logger::log_level = 4;
  ToplibSearchArgument argu_processor;
  if (!argu_processor.parse(argc, argv)) {
    return argu_processor.helpShown() ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  std::map<std::string, std::string> arguments = argu_processor.getArguments();
  ToplibSearchArgument::outputArguments(std::cout, arguments);
  auto start = std::chrono::steady_clock::now();

  std::cout << "Reading the library - started." << std::endl;
  toplib::LibSearchLibrary lib =
      toplib::lib_search::readLibrary(arguments["libraryFileName"]);
  int target_num = 0;
  for (double flag : lib.flags) {
    target_num += flag == 1 ? 1 : 0;
  }
  std::cout << "Reading the library - finished: " << lib.row_num
            << " identified spectra (" << target_num << " targets, "
            << lib.row_num - target_num << " decoys) of "
            << lib.table_row_num << "." << std::endl;

  std::cout << "Reading spectra - started." << std::endl;
  toplib::LibData queries =
      toplib::lib_msalign_reader::readMsalign(arguments["msalignFileName"]);
  std::cout << "Reading spectra - finished: " << queries.spectra.size()
            << " spectra with masses and a precursor mass." << std::endl;

  std::cout << "Searching the library - started." << std::endl;
  toplib::LibSearchParams params;
  params.tol_type = arguments["precErrorType"] == "Da"
                        ? toplib::LibTolType::kDa
                        : toplib::LibTolType::kPpm;
  params.precursor_tol =
      std::strtod(arguments["precErrorTolerance"].c_str(), nullptr);
  params.fragment_ppm =
      std::strtod(arguments["fragErrorTolerance"].c_str(), nullptr);
  params.charge_match = arguments["chargeMatch"] == "yes";
  toplib::LibSearchResult result =
      toplib::lib_search::search(lib, queries, params);
  int match_target_num = 0;
  for (const toplib::LibMatch& m : result.matches) {
    match_target_num += lib.flags[m.lib_row] == 1 ? 1 : 0;
  }
  std::cout << "Searching the library - finished: " << result.searched_num
            << " spectra with a precursor intensity and at least two masses, "
            << result.matches.size() << " with a matching library spectrum ("
            << match_target_num << " targets, "
            << result.matches.size() - match_target_num << " decoys)."
            << std::endl;

  if (!arguments["matchFileName"].empty()) {
    toplib::lib_search::writeMatches(arguments["matchFileName"], lib, queries,
                                     result);
  }

  if (result.cutoff_found) {
    std::cout << "FDR - cosine distance cutoff "
              << toplib::lib_numeric::pythonRepr(result.cutoff) << " (cutoff "
              << result.cutoff_index << "): " << result.cutoff_target_num
              << " targets, " << result.cutoff_decoy_num << " decoys, FDR "
              << toplib::lib_numeric::pythonRepr(result.cutoff_fdr) << "."
              << std::endl;
  } else {
    LOG_WARN("No match has a cosine distance below 0.999; no FDR cutoff.");
  }

  std::cout << "Writing the identifications - started." << std::endl;
  toplib::lib_search::writeResult(arguments["outputFileName"], lib, queries,
                                  result);
  std::cout << "Writing the identifications - finished: "
            << result.reported.size()
            << " target matches with a similarity of at least 0.3 in "
            << arguments["outputFileName"] << "." << std::endl;
  auto end = std::chrono::steady_clock::now();
  double secs = std::chrono::duration<double>(end - start).count();
  std::cout << "Running time: " << secs << " seconds." << std::endl;
  return EXIT_SUCCESS;
}
