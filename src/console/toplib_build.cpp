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
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "common/util/logger.hpp"
#include "console/toplib_build_argument.hpp"
#include "toplib/lib_builder.hpp"
#include "toplib/lib_msalign_reader.hpp"
#include "toplib/lib_spectrum.hpp"

using namespace toppic;

int main(int argc, char* argv[]) {
  toppic::logger::log_level = 4;
  ToplibBuildArgument argu_processor;
  if (!argu_processor.parse(argc, argv)) {
    return argu_processor.helpShown() ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  std::map<std::string, std::string> arguments = argu_processor.getArguments();
  ToplibBuildArgument::outputArguments(std::cout, arguments);
  auto start = std::chrono::steady_clock::now();

  std::cout << "Reading spectra - started." << std::endl;
  bool python_compatible = arguments["pythonCompatible"] == "true";
  toplib::LibData data;
  if (arguments["prsmFileName"].empty()) {
    data = toplib::lib_msalign_reader::readLibraryMsalign(
        arguments["msalignFileName"]);
  } else {
    data = toplib::lib_msalign_reader::readMsalignAndPrsms(
        arguments["msalignFileName"], arguments["prsmFileName"],
        python_compatible);
  }
  int identified = 0;
  for (const toplib::LibSpectrum& sp : data.spectra) {
    identified += sp.id_values.empty() ? 0 : 1;
  }
  std::cout << "Reading spectra - finished: " << data.spectra.size()
            << " spectra, " << identified << " identified." << std::endl;

  double cutoff = std::strtod(arguments["ppirCutoff"].c_str(), nullptr);
  if (cutoff != 0) {
    size_t before = data.spectra.size();
    toplib::lib_builder::filterPpir(data, cutoff);
    std::cout << "PPIR filter - " << data.spectra.size() << " of " << before
              << " spectra kept." << std::endl;
  }

  std::cout << "Clustering spectra and building representatives - started."
            << std::endl;
  toplib::LibRepMethod method = arguments["representative"] == "single"
                                    ? toplib::LibRepMethod::kSingle
                                    : toplib::LibRepMethod::kAverage;
  bool disable_different_proteoform_merge =
      arguments["disableDifferentProteoformMerge"] == "true";
  toplib::LibLibrary lib = toplib::lib_builder::build(
      data, method, disable_different_proteoform_merge, python_compatible);
  int written = 0;
  for (const toplib::LibRep& rep : lib.reps) {
    written += rep.identified ? 1 : 0;
  }
  std::cout << "Clustering spectra and building representatives - finished: "
            << lib.valid_spectrum_num << " spectra with at least two masses and "
            << "a precursor intensity, " << lib.precursor_cluster_num
            << " precursor clusters, " << lib.reps.size() << " clusters ("
            << lib.unidentified_cluster_num << " without an identification, "
            << lib.inconsistent_cluster_num
            << " losing it to a better cluster of the same proteoform and "
            << "charge)." << std::endl;

  if (!arguments["clusterFileName"].empty()) {
    std::ofstream out(arguments["clusterFileName"]);
    out << "cluster_id\tspectrum_index\tspectrum_id" << std::endl;
    for (size_t k = 0; k < lib.cluster_members.size(); k++) {
      for (int s : lib.cluster_members[k]) {
        out << k << "\t" << s << "\t" << data.spectra[s].spectrum_id
            << std::endl;
      }
    }
  }

  std::cout << "Writing the library - started." << std::endl;
  toplib::lib_sql_writer::write(lib, arguments["outputFileName"]);
  std::cout << "Writing the library - finished: " << written << " targets and "
            << written << " decoys in " << arguments["outputFileName"] << "."
            << std::endl;
  auto end = std::chrono::steady_clock::now();
  double secs = std::chrono::duration<double>(end - start).count();
  std::cout << "Running time: " << secs << " seconds." << std::endl;
  return EXIT_SUCCESS;
}
