// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iomanip>
#include <string>
#include <algorithm>

#include "boost/algorithm/string.hpp"

#include "base/file_util.hpp"
#include "base/xml_dom_util.hpp"
#include "feature/deconv_argu.hpp"

namespace prot {

DeconvArgument::DeconvArgument() {
  initArguments();
}

void DeconvArgument::initArguments() {
  arguments_["executiveDir"] = "";
  arguments_["spectrumFileName"] = "";
  arguments_["inputType"] = "mzXML";
  arguments_["outputType"] = "msalign";
  arguments_["refinePrecMass"]="true";
  arguments_["missingLevelOne"] = "false";
  arguments_["maxCharge"] = "30";
  arguments_["maxMass"] = "100000";
  arguments_["mzError"] = "0.02";
  arguments_["snRatio"] = "1.0";
  arguments_["keepUnusedPeaks"] = "false";
  arguments_["outMultipleMass"] = "false";
  arguments_["precWindow"] = "2.0";
}

void DeconvArgument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: toppfd [options] spectrum-file-name" << std::endl;
  std::cout << desc << std::endl;
}

bool DeconvArgument::parse(int argc, char* argv[]) {
  std::string spectrum_file_name = "";
  std::string output_type = "";
  std::string max_charge = "";
  std::string max_mass = "";
  std::string mz_error = "";
  std::string sn_ratio = "";

  // Define and parse the program options
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options()
        ("help,h", "Print this help message.")
        ("output,o",po::value<std::string>(&output_type),
         "<mgf|text|msalign>. Output file format: mgf, text or msalign. Default format is msalign.")
        ("max-charge,c", po::value<std::string> (&max_charge),
         "<integer value>. Set the maximum charge state of fragment and precursor ions. Default value is 30.")
        ("max-mass,m", po::value<std::string> (&max_mass),
         "<float value>. Set the maximum monoisopotic mass of fragment ions. Default value is 100000 Dalton.")
        ("mz-error,e", po::value<std::string> (&mz_error),
         "<float value>. Set the error tolerance of m/z values of peaks. Default value is 0.02 Dalton.")
        ("sn-ratio,s", po::value<std::string> (&sn_ratio),
         "<float value>. Set the signal/noise ratio. Default value is 1.")
        ("missing-level-one,n","The input spectrum file does not contain MS1 spectra.")
        ("multiple-mass,u", "Output multiple masses for one envelope.")
        ;
    po::options_description desc("Options");

    desc.add_options() 
        ("help,h", "Print this help message.") 
        ("output,o",po::value<std::string>(&output_type),
         "<mgf|text|msalign>. Output file format: mgf, text or msalign. Default format is msalign.")
        ("max-charge,c", po::value<std::string> (&max_charge),
         "<integer value>. Set the maximum charge state of fragment and precursor ions. Default value is 30.")
        ("max-mass,m", po::value<std::string> (&max_mass),
         "<float value>. Set the maximum monoisopotic mass of fragment ions. Default value is 100000 Dalton.")
        ("mz-error,e", po::value<std::string> (&mz_error),
         "<float value>. Set the error tolerance of m/z values of peaks. Default value is 0.02 Dalton.")
        ("sn-ratio,s", po::value<std::string> (&sn_ratio),
         "<float value>. Set the signal/noise ratio. Default value is 1.")
        ("missing-level-one,n","The input spectrum file does not contain MS1 spectra.")
        ("multiple-mass,u", "Output multiple masses for one envelope.")
        ("spectrum-file-name", po::value<std::string>(&spectrum_file_name)->required(), "Spectrum file name with its path.")
        ("keep,k", "Report monoisotopic masses extracted from low quality isotopic envelopes.")
        ;

    po::positional_options_description positional_options;
    positional_options.add("spectrum-file-name", 1);


    std::string app_name;
    //= boost::filesystem::basename(argv[0]);
    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(), vm);
      if ( vm.count("help") ) {
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
    std::string argv_0(argv[0]);
    arguments_["executiveDir"] = FileUtil::getExecutiveDir(argv_0);
    arguments_["spectrumFileName"] = spectrum_file_name;
    if (vm.count("output")) {
      arguments_["outputType"] = output_type;
    }
    if (vm.count("max-charge")) {
      arguments_["maxCharge"] = max_charge;
    }
    if (vm.count("keep")) {
      arguments_["keepUnusedPeaks"] = "true";
    }
    if (vm.count("max-mass")) {
      arguments_["maxMass"] = max_mass;
    }
    if (vm.count("mz-error")) {
      arguments_["mzError"] = mz_error;
    }
    if (vm.count("sn-ratio")) {
      arguments_["snRatio"] = sn_ratio;
    }
    if (vm.count("missing-level-one")) {
      arguments_["missingLevelOne"] = "true";
    }
    if (vm.count("multiple-mass")) {
      arguments_["outMultipleMass"] = "true";
    }
  }
  catch(std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line "
        << e.what() << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool DeconvArgument::validateArguments() {
  if (!boost::filesystem::exists(arguments_["spectrumFileName"])) {
    LOG_ERROR("Spectrum file " << arguments_["spectrumFileName"] << " does not exist!");
    return false;
  }
  std::string output_type = arguments_["outputType"];
  std::transform(output_type.begin(), output_type.end(), output_type.begin(), ::tolower);
  if (output_type != "msalign" && output_type != "mgf" && output_type != "text") {
    LOG_ERROR("Output type " << output_type << " error! The value should be mgf|text|msalign!");
    return false;
  }
  return true;
}

}  // namespace prot
