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

#include "console/toplib_search_argument.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "common/util/logger.hpp"
#include "common/util/version.hpp"

namespace toppic {

namespace {

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

bool isNonNegativeNumber(const std::string& s) {
  char* end = nullptr;
  double v = std::strtod(s.c_str(), &end);
  return end != s.c_str() && *end == '\0' && std::isfinite(v) && v >= 0;
}

}  // namespace

ToplibSearchArgument::ToplibSearchArgument() { arguments_ = initArguments(); }

std::map<std::string, std::string> ToplibSearchArgument::initArguments() {
  std::map<std::string, std::string> arguments;
  arguments["libraryFileName"] = "";
  arguments["msalignFileName"] = "";
  arguments["precErrorType"] = "ppm";
  // empty: 10 ppm or 2.2 Da, set in validateArguments
  arguments["precErrorTolerance"] = "";
  arguments["fragErrorTolerance"] = "10";
  arguments["chargeMatch"] = "no";
  arguments["outputFileName"] = "toplib_output/query_res.tsv";
  arguments["matchFileName"] = "";
  return arguments;
}

void ToplibSearchArgument::outputArguments(
    std::ostream& output, std::map<std::string, std::string> arguments) {
  int gap = 34;
  output << "******************** Parameters ********************" << std::endl;
  output << std::setw(gap) << std::left << "Library file:" << "\t"
         << arguments["libraryFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Spectrum file:" << "\t"
         << arguments["msalignFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Precursor error tolerance:" << "\t"
         << arguments["precErrorTolerance"] << " " << arguments["precErrorType"]
         << std::endl;
  output << std::setw(gap) << std::left << "Fragment error tolerance:" << "\t"
         << arguments["fragErrorTolerance"] << " ppm" << std::endl;
  output << std::setw(gap) << std::left << "Precursor charge matching:" << "\t"
         << arguments["chargeMatch"] << std::endl;
  output << std::setw(gap) << std::left << "Output file:" << "\t"
         << arguments["outputFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Version:" << "\t"
         << Version::getVersion() << std::endl;
  output << "******************** Parameters ********************" << std::endl;
}

void ToplibSearchArgument::showUsage(
    boost::program_options::options_description& desc) {
  std::cout << "Usage: toplib_search [options] <library file> <msalign file>"
            << std::endl;
  std::cout << desc << std::endl;
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool ToplibSearchArgument::parse(int argc, char* argv[]) {
  std::string prec_type;
  std::string prec_tol;
  std::string frag_tol;
  std::string charge;
  std::string output;
  std::string match_file;
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");
    display_desc.add_options()("help,h", "Print the help message.")(
        "precursor-error-type,T", po::value<std::string>(&prec_type),
        "<ppm|Da>. The unit of the precursor mass error tolerance. Default: "
        "ppm.")(
        "precursor-error,E", po::value<std::string>(&prec_tol),
        "<a nonnegative number>. The precursor mass error tolerance. In "
        "ppm (of the library spectrum's precursor mass), a library spectrum "
        "matches when its precursor mass is within the tolerance of the "
        "query's or one isotope shift (1.00235 Da) away from it within the "
        "tolerance; in Da, when it is within the tolerance. Default: 10 ppm "
        "or 2.2 Da.")(
        "fragment-error,e", po::value<std::string>(&frag_tol),
        "<a nonnegative number>. The fragment mass error tolerance in ppm. "
        "Default: 10.")(
        "charge-state,c", po::value<std::string>(&charge),
        "<yes|no>. Whether a library spectrum must have the query's "
        "precursor charge. Default: no.")(
        "output,o", po::value<std::string>(&output),
        "<file name>. The table of identifications. Default: "
        "toplib_output/query_res.tsv.");
    po::options_description desc("Options");
    desc.add(display_desc);
    // hidden: every query's best library spectrum before the FDR step, for
    // comparing with the Python script
    desc.add_options()("match-file", po::value<std::string>(&match_file), "")(
        "input-file", po::value<std::vector<std::string>>()->multitoken(),
        "The library file and the msalign file.");
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
    if (files.size() != 2) {
      std::cerr << "ERROR: give a library file and an msalign file."
                << std::endl
                << std::endl;
      showUsage(display_desc);
      return false;
    }
    arguments_["libraryFileName"] = files[0];
    arguments_["msalignFileName"] = files[1];
    if (vm.count("precursor-error-type")) {
      arguments_["precErrorType"] = prec_type;
    }
    if (vm.count("precursor-error")) {
      arguments_["precErrorTolerance"] = prec_tol;
    }
    if (vm.count("fragment-error")) {
      arguments_["fragErrorTolerance"] = frag_tol;
    }
    if (vm.count("charge-state")) {
      arguments_["chargeMatch"] = charge;
    }
    if (vm.count("output")) {
      arguments_["outputFileName"] = output;
    }
    if (vm.count("match-file")) {
      arguments_["matchFileName"] = match_file;
    }
  } catch (std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line" << e.what()
              << ", application will now exit" << std::endl;
    return false;
  }
  return validateArguments();
}

bool ToplibSearchArgument::validateArguments() {
  for (const char* key : {"libraryFileName", "msalignFileName"}) {
    const std::string& f = arguments_[key];
    if (!std::filesystem::exists(f)) {
      LOG_ERROR(f << " does not exist!");
      return false;
    }
  }
  std::string type = lower(arguments_["precErrorType"]);
  if (type != "ppm" && type != "da") {
    LOG_ERROR("The precursor error type must be ppm or Da, not "
              << arguments_["precErrorType"] << "!");
    return false;
  }
  arguments_["precErrorType"] = type == "da" ? "Da" : "ppm";
  if (arguments_["precErrorTolerance"].empty()) {
    arguments_["precErrorTolerance"] = type == "da" ? "2.2" : "10";
  }
  for (const char* key : {"precErrorTolerance", "fragErrorTolerance"}) {
    if (!isNonNegativeNumber(arguments_[key])) {
      LOG_ERROR("The error tolerance " << arguments_[key]
                                       << " is not a nonnegative number!");
      return false;
    }
  }
  std::string charge = lower(arguments_["chargeMatch"]);
  if (charge != "yes" && charge != "no") {
    LOG_ERROR("The charge state option must be yes or no, not "
              << arguments_["chargeMatch"] << "!");
    return false;
  }
  arguments_["chargeMatch"] = charge;
  return true;
}

}  // namespace toppic
