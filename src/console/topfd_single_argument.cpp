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

#include <iostream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "console/topfd_argument.hpp"

namespace toppic {

Argument::Argument() {
  initArguments();
}

void Argument::initArguments() {
  arguments_["executiveDir"] = "";
  arguments_["resourceDir"] = "";
  arguments_["precursorMass"] = "100000";
  arguments_["precursorCharge"] = "30";
  arguments_["mzError"] = "0.02";
  arguments_["msTwoSnRatio"] = "1.0";
  arguments_["keepUnusedPeaks"] = "false";
  arguments_["doFinalFiltering"] = "true";
  arguments_["outputMatchEnv"] = "false";
}

void Argument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: toppfd [options] spectrum-file-name" << std::endl;
  std::cout << desc << std::endl;
}

bool Argument::parse(int argc, char* argv[]) {
  std::string prec_charge = "";
  std::string prec_mass = "";
  std::string mz_error = "";
  std::string ms_two_sn_ratio = "";

  // Define and parse the program options
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options()
        ("help,h", "Print this help message.")
        ("prec-charge,c", po::value<std::string> (&prec_charge),
         "<a positive integer>. Set the precursor charge state. The default value is 30.")
        ("prec-mass,m", po::value<std::string> (&prec_mass),
         "<a positive number>. Set the precusor monoisotopic mass. The default value is 100000 Dalton.")
        ("mz-error,e", po::value<std::string> (&mz_error),
         "<a positive number>. Set the error tolerance of m/z values of spectral peaks. The default value is 0.02 m/z.")
        ("ms-two-sn-ratio,t", po::value<std::string> (&ms_two_sn_ratio),
         "<a positive number>. Set the signal/noise ratio for MS/MS spectra. The default value is 1.")
        ;
    po::options_description desc("Options");

    desc.add_options() 
        ("help,h", "Print this help message.") 
        ("prec-charge,c", po::value<std::string> (&prec_charge), "")
        ("prec-mass,m", po::value<std::string> (&prec_mass), "")
        ("mz-error,e", po::value<std::string> (&mz_error), "")
        ("ms-two-sn-ratio,t", po::value<std::string> (&ms_two_sn_ratio), "")
        ("keep,k", "Report monoisotopic masses extracted from low quality isotopic envelopes.")
        ("spectrum-file-name", po::value<std::vector<std::string> >()->multitoken()->required(), "Spectrum file name with its path.")
        ;

    po::positional_options_description positional_options;
    positional_options.add("spectrum-file-name", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(), vm);
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

    // get the execution directory
    std::string argv_0(argv[0]);
    arguments_["executiveDir"] = file_util::getExecutiveDir(argv_0);
    if (file_util::checkSpace(arguments_["executiveDir"])) {
      LOG_ERROR("Current directory " << arguments_["executiveDir"] << " contains space and will cause errors in the program!")
      exit(EXIT_FAILURE);
    }

    arguments_["resourceDir"] = file_util::getResourceDir(arguments_["executiveDir"]);

    if (vm.count("prec-charge")) {
      arguments_["precursorCharge"] = prec_charge;
    }

    if (vm.count("keep")) {
      arguments_["keepUnusedPeaks"] = "true";
    }

    if (vm.count("prec-mass")) {
      arguments_["precursorMass"] = prec_mass;
    }

    if (vm.count("mz-error")) {
      arguments_["mzError"] = mz_error;
    }

    if (vm.count("ms-two-sn-ratio")) {
      arguments_["msTwoSnRatio"] = ms_two_sn_ratio;
    }

    if (vm.count("spectrum-file-name")) {
      spec_file_list_ = vm["spectrum-file-name"].as<std::vector<std::string> >(); 
    }

  }
  catch(std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line "
        << e.what() << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool Argument::validateArguments() {
  if (!file_util::exists(arguments_["resourceDir"])) {
    LOG_ERROR("Resource direcotry " << arguments_["resourceDir"] << " does not exist!\nPlease check if file directory contains "
    << "unproper characters such as spaces/quotation makrks");
    return false;
  }

  for (size_t k = 0; k < spec_file_list_.size(); k++) {
    if (!file_util::exists(spec_file_list_[k])) {
      LOG_ERROR(spec_file_list_[k] << " does not exist!\nPlease check if file directory contains unproper characters such as spaces/quotation makrks");
      return false;
    }
  }
  return true;
}

}  // namespace toppic
