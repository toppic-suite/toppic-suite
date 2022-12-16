//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <iomanip>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "console/topdiff_argument.hpp"
#include "common/util/version.hpp"

namespace toppic {

TopDiffArgument::TopDiffArgument() {
  arguments_ = initArguments();
}

std::map<std::string, std::string> TopDiffArgument::initArguments() {
  std::map<std::string, std::string> arguments;
  arguments["errorTolerance"] = "1.2";
  arguments["toolName"] = "toppic";
  arguments["mergedOutputFileName"] = "sample_diff.tsv";
  arguments["version"] = "";
  return arguments;
}

void TopDiffArgument::outputArguments(std::ostream &output, 
                                      const std::string &sep,
                                      std::map<std::string, std::string> arguments) {
  int gap = 26;
  output << "******************** Parameters ********************" << std::endl;
  output << std::setw(gap) << std::left << "Error tolerance:" << sep << arguments["errorTolerance"] << " Dalton " << std::endl;
  output << std::setw(gap) << std::left << "Database search tool name:" << sep << arguments["toolName"] << std::endl;
  output << std::setw(gap) << std::left << "Output file name:" << sep << arguments["mergedOutputFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Version:" << sep << arguments["version"] << std::endl;
  output << "******************** Parameters ********************" << std::endl;
}

void TopDiffArgument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: topdiff [options] spectrum-file-names" << std::endl; 
  std::cout << desc << std::endl; 
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool TopDiffArgument::parse(int argc, char* argv[]) {
  std::string error_tole = "1.2";
  std::string tool_name = "toppic";
  std::string merged_output_name = "sample_diff.tsv";

  /** Define and parse the program options*/
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");
    display_desc.add_options() 
        ("help,h", "Print the help message.")
        ("error-tolerance,e", po::value<std::string>(&error_tole) , "Error tolerance of precursor masses for mapping identified proteoforms. Default: 1.2 Dalton.")
        ("tool-name,t", po::value<std::string>(&tool_name) , "<toppic|topmg>. Database search tool name: toppic or topmg. Default: toppic.")
        ("output,o", po::value<std::string>(&merged_output_name) , "Output file name. Default: sample_diff.tsv.");
    po::options_description desc("Options");

    desc.add_options() 
        ("help,h", "") 
        ("error-tolerance,e", po::value<std::string>(&error_tole) , "")
        ("tool-name,t", po::value<std::string>(&tool_name), "")
        ("output,o", po::value<std::string>(&merged_output_name) , "")
        ("spectrum-file-name", po::value<std::vector<std::string> >()->multitoken()->required(), "Spectrum file name with its path.");

    po::positional_options_description positional_options;
    positional_options.add("spectrum-file-name", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(),vm); 
      if (vm.count("help")) {
        showUsage(display_desc);
        return false;
      }
      po::notify(vm); 
      // throws on error, so do after help in case there are any problems
    } 
    catch(boost::program_options::required_option& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }
    catch(boost::program_options::error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }
    std::string argv_0 (argv[0]);
    arguments_["executiveDir"] = file_util::getExecutiveDir(argv_0);
    if (file_util::checkSpace(arguments_["executiveDir"])) {
      LOG_ERROR("Current directory " << arguments_["executiveDir"] << " contains space and will cause errors in the program!")
      exit(EXIT_FAILURE);
    }
    LOG_DEBUG("Executive Dir " << arguments_["executiveDir"]);
    arguments_["resourceDir"] = file_util::getResourceDir(arguments_["executiveDir"]);

    if (vm.count("spectrum-file-name")) {
      spectrum_file_list_ = vm["spectrum-file-name"].as<std::vector<std::string> >(); 
    }

    if (vm.count("error-tolerance")) {
      arguments_["errorTolerance"] = error_tole;
    }

    if (vm.count("tool-name")) {
      arguments_["toolName"] = tool_name;
    }

    if (vm.count("output")) {
      arguments_["mergedOutputFileName"] = merged_output_name; 
    }
  }
  catch(std::exception&e ) {
    std::cerr << "Unhandled Exception in parsing command line" << e.what() << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool TopDiffArgument::validateArguments() {
  if (!file_util::exists(arguments_["resourceDir"])) {
    LOG_ERROR("Resource direcotry " << arguments_["resourceDir"] << " does not exist!\n" 
              << "Please check if the file directory or name contains special characters such as spaces or quotation marks.");
    return false;
  }

  for (size_t k = 0; k < spectrum_file_list_.size(); k++) {
    if (!file_util::exists(spectrum_file_list_[k])) {
      LOG_ERROR(spectrum_file_list_[k] << " does not exist!\n"
                << "Please check if the file directory or name contains special characters such as spaces or quotation marks.");
      return false;
    }
  }
  if (spectrum_file_list_.size() < 2) {
    LOG_ERROR("Only 1 spectrum file added. Please add at least 2 spectum files to make a comparison.");
    return false;
  }

  return true;
}

}  // namespace toppic
