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

#include "console/topindex_argument.hpp"

#include <exception>
#include <filesystem>
#include <iomanip>
#include <string>
#include <utility>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/str_util.hpp"
#include "common/util/version.hpp"
#include "seq/fasta_util.hpp"

namespace toppic {

TopIndexArgument::TopIndexArgument() { arguments_ = initArguments(); }

std::map<std::string, std::string> TopIndexArgument::initArguments() {
  std::map<std::string, std::string> arguments;
  arguments["oriDatabaseFileName"] = "";
  arguments["databaseFileName"] = "";
  arguments["databaseBlockSize"] = "60000000";
  arguments["maxFragmentLength"] = "500";
  arguments["minBlockNum"] = "10";
  arguments["searchType"] = "TARGET";
  arguments["fixedMod"] = "";
  arguments["precErrorTolerance"] = "10";
  arguments["fragErrorTolerance"] = "10";
  arguments["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments["executiveDir"] = ".";
  arguments["resourceDir"] = "";
  arguments["threadNumber"] = "1";
  arguments["indexFormat"] = "2";
  arguments["version"] = "";

  // the following two arguments are used for the initializatio of prsm para
  arguments["groupSpectrumNumber"] = "1";
  arguments["activation"] = "FILE";
  return arguments;
}

void TopIndexArgument::outputArguments(
    std::ostream& output, const std::string& sep,
    std::map<std::string, std::string> arguments) {
  int gap = 36;
  output << "********************** Parameters **********************"
         << std::endl;
  output << std::setw(gap) << std::left << "Protein database file:" << sep
         << arguments["oriDatabaseFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Search type: " << sep
         << arguments["searchType"] << std::endl;

  if (arguments["fixedMod"] == "") {
    output << std::setw(gap) << std::left << "Fixed modifications:" << sep
           << "None" << std::endl;
  } else if (arguments["fixedMod"] == "C57") {
    output << std::setw(gap) << std::left << "Fixed modifications:" << sep
           << "C57:carbamidomethylation on cysteine" << std::endl;
  } else if (arguments["fixedMod"] == "C58") {
    output << std::setw(gap) << std::left << "Fixed modifications:" << sep
           << "C58:carboxymethylation on cysteine" << std::endl;
  } else {
    output << std::setw(gap) << std::left << "Fixed modifications:" << sep
           << arguments["fixedMod"] << std::endl;
  }

  output << std::setw(gap) << std::left
         << "Error tolerance for precursor masses:" << sep
         << arguments["precErrorTolerance"] << " ppm" << std::endl;
  output << std::setw(gap) << std::left
         << "Error tolerance for fragment masses:" << sep
         << arguments["fragErrorTolerance"] << " ppm" << std::endl;

  output << std::setw(gap) << std::left << "Allowed N-terminal forms:" << sep
         << arguments["allowProtMod"] << std::endl;
  output << std::setw(gap) << std::left << "Thread number:" << sep
         << arguments["threadNumber"] << std::endl;
  output << std::setw(gap) << std::left << "Version:" << sep
         << arguments["version"] << std::endl;
  output << "********************** Parameters **********************"
         << std::endl;
}

void TopIndexArgument::showUsage(
    boost::program_options::options_description& desc) {
  std::cout << "Usage: topindex [options] database-file-name" << std::endl;
  std::cout << desc << std::endl;
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool TopIndexArgument::parse(int argc, char* argv[]) {
  std::string database_file_name = "";
  std::string argument_file_name = "";
  std::string fixed_mod = "";
  std::string allow_mod = "";
  std::string frag_error_tole = "";
  std::string prec_error_tole = "";
  std::string thread_number = "";
  std::string index_format = "";

  /** Define and parse the program options*/
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options()("help,h", "Print the help message.")(
        "fixed-mod,f", po::value<std::string>(&fixed_mod),
        "<C57|C58|a fixed modification file>. Fixed modifications. Three "
        "available options: C57, C58, or the name of a text file containing "
        "the information of fixed modifications. When C57 is selected, "
        "carbamidomethylation on cysteine is the only fixed modification. When "
        "C58 is selected, carboxymethylation on cysteine is the only fixed "
        "modification.")(
        "n-terminal-form,n", po::value<std::string>(&allow_mod),
        "<a list of allowed N-terminal forms>. N-terminal forms of proteins. "
        "Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, "
        "and M_ACETYLATION. NONE stands for no modifications, NME for "
        "N-terminal methionine excision, NME_ACETYLATION for N-terminal "
        "acetylation after the initiator methionine is removed, and "
        "M_ACETYLATION for N-terminal methionine acetylation. When multiple "
        "forms are allowed, they are separated by commas. Default value: "
        "NONE,NME,NME_ACETYLATION,M_ACETYLATION.")(
        "decoy,d",
        "Index a database of the target proteins and shuffled decoy "
        "proteins, for TopPIC and TopMG searches with -d (--decoy).")(
        "frag-error-tolerance,e", po::value<std::string>(&frag_error_tole),
        "<a positive integer>. Error tolerance for fragment masses in PPM. "
        "The index files do not depend on it, but their names include it: "
        "give the -e value of the TopPIC or TopMG searches that will use "
        "them. Default value: 10.")(
        "prec-error-tolerance,E", po::value<std::string>(&prec_error_tole),
        "<a positive integer>. Error tolerance for precursor masses in PPM. "
        "Neither the index files nor their names depend on it. Default "
        "value: 10.")(
        "thread-number,u", po::value<std::string>(&thread_number),
        "<a positive integer>. Number of threads used in the computation. "
        "Default value: 1.");

    po::options_description desc("Options");

    desc.add_options()("help,h", "")("fixed-mod,f",
                                     po::value<std::string>(&fixed_mod), "")(
        "n-terminal-form,n", po::value<std::string>(&allow_mod), "")(
        "decoy,d", "")("frag-error-tolerance,e",
                       po::value<std::string>(&frag_error_tole), "")(
        "prec-error-tolerance,E", po::value<std::string>(&prec_error_tole), "")(
        "thread-number,u", po::value<std::string>(&thread_number), "")(
        "index-format", po::value<std::string>(&index_format), "")(
        "database-file-name",
        po::value<std::string>(&database_file_name)->required(),
        "Database file name with its path.");

    po::positional_options_description positional_options;
    positional_options.add("database-file-name", 1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv)
                    .options(desc)
                    .positional(positional_options)
                    .run(),
                vm);
      if (vm.count("help")) {
        showUsage(display_desc);
        return false;
      }
      po::notify(vm);
      // throws on error, so do after help in case there are any problems
    } catch (boost::program_options::required_option& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    } catch (boost::program_options::error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }
    std::string argv_0(argv[0]);
    if (vm.count("full-binary-path")) {
      arguments_["executiveDir"] = argv[0];
    } else {
      arguments_["executiveDir"] = file_util::getExecutiveDir(argv_0);
      if (file_util::checkSpace(arguments_["executiveDir"])) {
        LOG_ERROR("Current directory "
                  << arguments_["executiveDir"]
                  << " contains space and will cause errors in the program!");
        exit(EXIT_FAILURE);
      }
    }
    LOG_DEBUG("Executive Dir " << arguments_["executiveDir"]);

    arguments_["resourceDir"] =
        file_util::getResourceDir(arguments_["executiveDir"]);

    arguments_["oriDatabaseFileName"] = database_file_name;

    if (vm.count("decoy")) {
      arguments_["searchType"] = "TARGET+DECOY";
    }


    if (vm.count("fixed-mod")) {
      arguments_["fixedMod"] = fixed_mod;
    }

    if (vm.count("n-terminal-form")) {
      arguments_["allowProtMod"] = allow_mod;
    }

    if (vm.count("frag-error-tolerance")) {
      arguments_["fragErrorTolerance"] = frag_error_tole;
    }

    if (vm.count("prec-error-tolerance")) {
      arguments_["precErrorTolerance"] = prec_error_tole;
    }

    if (vm.count("thread-number")) {
      arguments_["threadNumber"] = thread_number;
    }
    // hidden: the index file format, 2 (zstd frames, the column bounds once;
    // the default) or 1 (the Boost binary archive of the earlier versions);
    // toppic and topmg read both (toppic_test experiment 150)
    if (vm.count("index-format")) {
      if (index_format != "1" && index_format != "2") {
        LOG_ERROR("The index format must be 1 or 2: " << index_format);
        return false;
      }
      arguments_["indexFormat"] = index_format;
    }
  } catch (std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line" << e.what()
              << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool TopIndexArgument::validateArguments() {
  if (!std::filesystem::exists(arguments_["resourceDir"])) {
    LOG_ERROR("Resource directory "
              << arguments_["resourceDir"] << " does not exist!\n"
              << "Please check if file directory or name contains special "
                 "characters such as spaces or quotation marks.");
    return false;
  }

  if (!std::filesystem::exists(arguments_["oriDatabaseFileName"])) {
    LOG_ERROR("Database file "
              << arguments_["oriDatabaseFileName"] << " does not exist!\n"
              << " Please check if file directory or name contains special "
                 "characters such as spaces or quotation marks.");
    return false;
  }

  if (!str_util::endsWith(arguments_["oriDatabaseFileName"], ".fasta") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".fa") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".FASTA") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".FA")) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"]
                               << " is not a fasta file!");
    return false;
  }

  if (arguments_["oriDatabaseFileName"].length() > 200) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"]
                               << " path is too long!");
    return false;
  }

  // A database whose protein ids (the header text up to the first
  // whitespace) repeat cannot be indexed by htslib: search its
  // _dedup.fasta copy, which keeps the first record of every id, instead.
  arguments_["oriDatabaseFileName"] =
      fasta_util::removeDuplicatedIds(arguments_["oriDatabaseFileName"]);

  if (arguments_["searchType"] == "TARGET+DECOY") {
    arguments_["databaseFileName"] =
        arguments_["oriDatabaseFileName"] + "_target_decoy";
  } else {
    arguments_["databaseFileName"] =
        arguments_["oriDatabaseFileName"] + "_target";
  }

  std::string search_type = arguments_["searchType"];
  if (search_type != "TARGET" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Search type "
              << search_type
              << " error! The value should be TARGET|TARGET+DECOY!");
    return false;
  }

  std::string allow_mod = arguments_["allowProtMod"];
  std::vector<std::string> strs = str_util::split(allow_mod, ",");
  for (size_t i = 0; i < strs.size(); i++) {
    if (strs[i] != "NONE" && strs[i] != "M_ACETYLATION" && strs[i] != "NME" &&
        strs[i] != "NME_ACETYLATION") {
      LOG_ERROR(
          "N-Terminal Variable PTM can only be NONE, M_ACETYLATION, NME or "
          "NME_ACETYLATION.");
      return false;
    }
  }

  for (const auto& [key, name] :
       {std::pair<const char*, const char*>{"precErrorTolerance",
                                            "Precursor error tolerance"},
        std::pair<const char*, const char*>{"fragErrorTolerance",
                                            "Fragment error tolerance"}}) {
    std::string tole_value = arguments_[key];
    try {
      double tole = std::stoi(tole_value);
      if (tole <= 0) {
        LOG_ERROR(name << ": " << tole_value
                       << " error! The value should be positive.");
        return false;
      }
    } catch (const std::exception& ex) {
      LOG_ERROR(name << ": " << tole_value << " should be a number.");
      return false;
    }
  }

  std::string thread_number = arguments_["threadNumber"];
  try {
    int num = std::stoi(thread_number.c_str());
    bool valid = mem_check::checkThreadNum(num, "topindex");
    if (!valid) {
      return false;
    }
  } catch (std::exception& e) {
    LOG_ERROR("Thread number " << thread_number << " should be a number.");
    return false;
  }

  return true;
}

}  // namespace toppic
