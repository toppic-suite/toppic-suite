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


#include "console/toplib_build_argument.hpp"

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>

#include "common/util/logger.hpp"
#include "common/util/version.hpp"

namespace toppic {

ToplibBuildArgument::ToplibBuildArgument() { arguments_ = initArguments(); }

std::map<std::string, std::string> ToplibBuildArgument::initArguments() {
  std::map<std::string, std::string> arguments;
  arguments["msalignFileName"] = "";
  arguments["prsmFileName"] = "";
  arguments["representative"] = "average";
  arguments["ppirCutoff"] = "0";
  arguments["outputFileName"] = "TopLib/lib_spectra_ms2.db";
  arguments["clusterFileName"] = "";
  arguments["pythonCompatible"] = "false";
  return arguments;
}

void ToplibBuildArgument::outputArguments(
    std::ostream& output, std::map<std::string, std::string> arguments) {
  int gap = 26;
  output << "******************** Parameters ********************" << std::endl;
  output << std::setw(gap) << std::left << "Spectrum file:" << "\t"
         << arguments["msalignFileName"] << std::endl;
  output << std::setw(gap) << std::left << "PrSM file:" << "\t"
         << (arguments["prsmFileName"].empty()
                 ? "none (identifications in the spectrum headers)"
                 : arguments["prsmFileName"])
         << std::endl;
  output << std::setw(gap) << std::left << "Representative:" << "\t"
         << arguments["representative"] << std::endl;
  output << std::setw(gap) << std::left << "PPIR cutoff:" << "\t"
         << arguments["ppirCutoff"] << std::endl;
  output << std::setw(gap) << std::left << "Output file:" << "\t"
         << arguments["outputFileName"] << std::endl;
  if (arguments["pythonCompatible"] == "true") {
    output << std::setw(gap) << std::left << "Python compatible:" << "\t"
           << "true" << std::endl;
  }
  output << std::setw(gap) << std::left << "Version:" << "\t"
         << Version::getVersion() << std::endl;
  output << "******************** Parameters ********************" << std::endl;
}

void ToplibBuildArgument::showUsage(
    boost::program_options::options_description& desc) {
  std::cout << "Usage: toplib_build [options] <msalign file> [<PrSM file>]"
            << std::endl;
  std::cout << desc << std::endl;
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool ToplibBuildArgument::parse(int argc, char* argv[]) {
  std::string representative;
  std::string ppir_cutoff;
  std::string output;
  std::string cluster_file;
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");
    display_desc.add_options()("help,h", "Print the help message.")(
        "representative,r", po::value<std::string>(&representative),
        "<average|single>. The representative spectrum of a cluster: the "
        "merged spectrum of its members (average) or its member with the "
        "lowest E-value (single). Default: average.")(
        "ppir-cutoff,c", po::value<std::string>(&ppir_cutoff),
        "<a number between 0 and 1>. Keep only the spectra whose first "
        "precursor has at least this fraction of the summed intensity of the "
        "spectrum's precursors (the precursor-to-precursor intensity ratio). "
        "Default: 0 (keep all spectra).")(
        "output,o", po::value<std::string>(&output),
        "<file name>. The library database. Default: "
        "TopLib/lib_spectra_ms2.db.");
    po::options_description desc("Options");
    desc.add(display_desc);
    // hidden: the cluster of every spectrum, for testing; the Python
    // script's behaviour (spectra of different proteoforms may share a
    // cluster, TopLib's PrSM table column names only), for comparing with it
    desc.add_options()("cluster-file", po::value<std::string>(&cluster_file),
                       "")("python-compatible", "")("input-file",
                           po::value<std::vector<std::string>>()->multitoken(),
                           "The msalign file and the PrSM file.");
    po::positional_options_description positional_options;
    positional_options.add("input-file", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv)
                    .options(desc)
                    .positional(positional_options)
                    .run(),
                vm);
      if (vm.count("help")) {
        showUsage(display_desc);
        help_shown_ = true;
        return false;
      }
      po::notify(vm);
    } catch (boost::program_options::error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }
    std::vector<std::string> files;
    if (vm.count("input-file")) {
      files = vm["input-file"].as<std::vector<std::string>>();
    }
    if (files.empty() || files.size() > 2) {
      std::cerr << "ERROR: give an msalign file and, optionally, a PrSM file."
                << std::endl
                << std::endl;
      showUsage(display_desc);
      return false;
    }
    arguments_["msalignFileName"] = files[0];
    if (files.size() == 2) {
      arguments_["prsmFileName"] = files[1];
    }
    if (vm.count("representative")) {
      arguments_["representative"] = representative;
    }
    if (vm.count("ppir-cutoff")) {
      arguments_["ppirCutoff"] = ppir_cutoff;
    }
    if (vm.count("output")) {
      arguments_["outputFileName"] = output;
    }
    if (vm.count("cluster-file")) {
      arguments_["clusterFileName"] = cluster_file;
    }
    if (vm.count("python-compatible")) {
      arguments_["pythonCompatible"] = "true";
    }
  } catch (std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line" << e.what()
              << ", application will now exit" << std::endl;
    return false;
  }
  return validateArguments();
}

bool ToplibBuildArgument::validateArguments() {
  for (const char* key : {"msalignFileName", "prsmFileName"}) {
    const std::string& f = arguments_[key];
    if (!f.empty() && !std::filesystem::exists(f)) {
      LOG_ERROR(f << " does not exist!");
      return false;
    }
  }
  const std::string& rep = arguments_["representative"];
  if (rep != "average" && rep != "single") {
    LOG_ERROR("The representative must be average or single, not " << rep
                                                                   << "!");
    return false;
  }
  const std::string& cutoff = arguments_["ppirCutoff"];
  char* end = nullptr;
  double v = std::strtod(cutoff.c_str(), &end);
  if (end == cutoff.c_str() || *end != '\0' || v < 0 || v > 1) {
    LOG_ERROR("The PPIR cutoff " << cutoff
                                 << " is not a number between 0 and 1!");
    return false;
  }
  return true;
}

}  // namespace toppic
