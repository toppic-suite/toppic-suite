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

#ifndef BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_SYSTEM_NO_DEPRECATED 1
#endif

#include <iomanip>

#include "boost/thread/thread.hpp"

#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"

#include "console/topindex_argument.hpp"

namespace toppic {

Argument::Argument() {
  initArguments();
}

void Argument::initArguments() {
  arguments_["oriDatabaseFileName"]="";
  arguments_["databaseFileName"] = "";
  arguments_["databaseBlockSize"] = "1000000";
  arguments_["searchType"] = "TARGET";
  arguments_["fixedMod"] = "";
  arguments_["massErrorTolerance"] = "15";
  arguments_["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments_["executiveDir"] = ".";
  arguments_["resourceDir"] = "";
  arguments_["threadNumber"] = "1";

  // filtering result number is for diagonal filter
  arguments_["filteringResultNumber"] = "20";
  
  // the following two arguments are used for the initializatio of prsm para
  arguments_["groupSpectrumNumber"] = "1";
  arguments_["activation"] = "FILE";
}

void Argument::outputArguments(std::ostream &output, 
                               std::map<std::string, std::string> arguments) {
  output << "********************** Parameters **********************" << std::endl;
  output << std::setw(44) << std::left << "Protein database file: " << "\t" << arguments["oriDatabaseFileName"] << std::endl;
  output << std::setw(44) << std::left << "Search type: " << "\t" << arguments["searchType"] << std::endl;

  if (arguments["fixedMod"] == "") {
    output << std::setw(44) << std::left << "Fixed modifications: " << "\t" << "None" << std::endl;
  } 
  else if (arguments["fixedMod"] == "C57") {
    output << std::setw(44) << std::left << "Fixed modifications: " << "\t" << "C57:carbamidomethylation on cysteine" << std::endl;
  }
  else if (arguments["fixedMod"] == "C58") {
    output << std::setw(44) << std::left << "Fixed modifications: " << "\t" << "C58:carboxymethylation on cysteine" << std::endl;
  }
  else {
    output << std::setw(44) << std::left << "Fixed modifications:," << arguments["fixedMod"] << std::endl;
  }

  output << std::setw(44) << std::left << "Maximum number of unexpected modifications: " << "\t" << arguments["ptmNumber"] << std::endl;
  output << std::setw(44) << std::left << "Error tolerance for matching masses: " << "\t" << arguments["massErrorTolerance"] << " ppm" << std::endl;

  output << std::setw(44) << std::left << "Allowed N-terminal forms: " << "\t" <<  arguments["allowProtMod"] << std::endl;
  output << std::setw(44) << std::left << "Thread number: " << "\t" << arguments["threadNumber"] << std::endl;
  output << "********************** Parameters **********************" << std::endl;
}


void Argument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: topindex [options] database-file-name" << std::endl; 
  std::cout << desc << std::endl; 
}

bool Argument::parse(int argc, char* argv[]) {
  
  std::string database_file_name = "";
  std::string argument_file_name = "";
  std::string fixed_mod = "";
  std::string allow_mod = "";
  std::string mass_error_tole = "";
  std::string thread_number = "";

  /** Define and parse the program options*/
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options() 
        ("help,h", "Print the help message.") 
        ("fixed-mod,f", po::value<std::string> (&fixed_mod), 
         "<C57|C58|a fixed modification file>. Fixed modifications. Three available options: C57, C58, or the name of a text file containing the information of fixed modifications. When C57 is selected, carbamidomethylation on cysteine is the only fixed modification. When C58 is selected, carboxymethylation on cysteine is the only fixed modification.")
        ("n-terminal-form,n", po::value<std::string> (&allow_mod), 
         "<a list of allowed N-terminal forms>. N-terminal forms of proteins. Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, and M_ACETYLATION. NONE stands for no modifications, NME for N-terminal methionine excision, NME_ACETYLATION for N-terminal acetylation after the initiator methionine is removed, and M_ACETYLATION for N-terminal methionine acetylation. When multiple forms are allowed, they are separated by commas. Default value: NONE,NME,NME_ACETYLATION,M_ACETYLATION.")
        ("decoy,d", "Use a decoy protein database to estimate false discovery rates.")
        ("mass-error-tolerance,e", po::value<std::string> (&mass_error_tole), "<a positive integer>. Error tolerance for precursor and fragment masses in PPM. Default value: 15.")
        ("thread-number,u", po::value<std::string> (&thread_number), "<a positive integer>. Number of threads used in the computation. Default value: 1.");

    po::options_description desc("Options");

    desc.add_options() 
        ("help,h", "") 
        ("fixed-mod,f", po::value<std::string> (&fixed_mod), "")
        ("n-terminal-form,n", po::value<std::string> (&allow_mod), "")
        ("decoy,d", "")
        ("mass-error-tolerance,e", po::value<std::string> (&mass_error_tole), "")
        ("thread-number,u", po::value<std::string> (&thread_number), "")
        ("database-file-name", po::value<std::string>(&database_file_name)->required(), "Database file name with its path.");

    po::positional_options_description positional_options;
    positional_options.add("database-file-name", 1);

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
    if (vm.count("full-binary-path")) {
      arguments_["executiveDir"] = argv[0];
    } else {
      arguments_["executiveDir"] = file_util::getExecutiveDir(argv_0);
    }
    LOG_DEBUG("Executive Dir " << arguments_["executiveDir"]);

    arguments_["resourceDir"] = file_util::getResourceDir(arguments_["executiveDir"]);

    arguments_["oriDatabaseFileName"] = database_file_name;

    if (vm.count("decoy")) {
      arguments_["searchType"] = "TARGET+DECOY";
    }

    if (arguments_["searchType"] == "TARGET+DECOY") {
      arguments_["databaseFileName"]=arguments_["oriDatabaseFileName"] + "_target_decoy";
    } else {
      arguments_["databaseFileName"]=arguments_["oriDatabaseFileName"] + "_target";
    }

    if (vm.count("fixed-mod")) {
      arguments_["fixedMod"] = fixed_mod;
    }

    if (vm.count("n-terminal-form")) {
      arguments_["allowProtMod"] = allow_mod;
    }    

    if (vm.count("mass-error-tolerance")) {
      arguments_["massErrorTolerance"] = mass_error_tole;
    }

    if (vm.count("thread-number")) {
      arguments_["threadNumber"] = thread_number;
    }
  }
  catch(std::exception&e ) {
    std::cerr << "Unhandled Exception in parsing command line" << e.what() << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool Argument::validateArguments() {
  if (!file_util::exists(arguments_["resourceDir"])) {
    LOG_ERROR("Resource direcotry " << arguments_["resourceDir"] << " does not exist!");
    return false;
  }

  if (!file_util::exists(arguments_["oriDatabaseFileName"])) {
    LOG_ERROR("Database file " << arguments_["databaseFileName"] << " does not exist!");
    return false;
  }

  if (!str_util::endsWith(arguments_["oriDatabaseFileName"], ".fasta") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".fa")) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"] << " is not a fasta file!");
    return false;
  }

  if (arguments_["oriDatabaseFileName"].length() > 200) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"] << " path is too long!");
    return false;
  }

  std::string search_type = arguments_["searchType"];
  if(search_type != "TARGET" && search_type != "TARGET+DECOY"){
    LOG_ERROR("Search type " << search_type << " error! The value should be TARGET|TARGET+DECOY!");
    return false;
  }

  std::string allow_mod = arguments_["allowProtMod"]; 
  std::vector<std::string> strs = str_util::split(allow_mod, ",");
  for (size_t i = 0; i < strs.size(); i++) {
    if (strs[i] != "NONE" && strs[i] != "M_ACETYLATION" && strs[i] != "NME" && strs[i] != "NME_ACETYLATION") {
      LOG_ERROR("N-Terminal Variable PTM can only be NONE, M_ACETYLATION, NME or NME_ACETYLATION.");
      return false;
    }
  }

  std::string mass_error_tole_value = arguments_["massErrorTolerance"];
  try {
    double tole = std::stoi(mass_error_tole_value);
    if (tole <= 0) {
      LOG_ERROR("Mass error tolerance: " << mass_error_tole_value << " error! The value should be positive.");
      return false;
    }
  } catch (int e) {
    LOG_ERROR("Mass error tolerance: " << mass_error_tole_value << " should be a number.");
    return false;
  }


  std::string thread_number = arguments_["threadNumber"];
  try {
    int num = std::stoi(thread_number.c_str());
    if (num <= 0) {
      LOG_ERROR("Thread number " << thread_number << " error! The value should be positive.");
      return false;
    }
    int n = static_cast<int>(boost::thread::hardware_concurrency());
    if(num > n){
      LOG_ERROR("Thread number " << thread_number << " error! The value is too large. Only " << n << " threads are supported.");
      return false;
    }
  } catch (int e) {
    LOG_ERROR("Thread number " << thread_number << " should be a number.");
    return false;
  }

  return true;
}

}  // namespace toppic
