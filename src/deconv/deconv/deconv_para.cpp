//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <map>
#include <string>
#include <sstream>
#include <iomanip>

#include "common/util/version.hpp"
#include "common/util/str_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/time_util.hpp"
#include "deconv/deconv/deconv_para.hpp"

namespace toppic {

DeconvPara::DeconvPara(std::map<std::string, std::string> &arguments) {
  data_file_name_ = arguments["spectrumFileName"];

  fraction_id_ = std::stoi(arguments["fractionId"]);
  
  resource_dir_ = arguments["resourceDir"];

  missing_level_one_ = (arguments["missingLevelOne"] == "true");
  
  max_charge_ = std::stoi(arguments["maxCharge"]);
  
  max_mass_ = std::stod(arguments["maxMass"]);
  
  tolerance_ = std::stod(arguments["mzError"]);
  
  ms_two_sn_ratio_ = std::stod(arguments["msTwoSnRatio"]);
  
  ms_one_sn_ratio_ = std::stod(arguments["msOneSnRatio"]);
  
  keep_unused_peaks_ = (arguments["keepUnusedPeaks"] == "true");
  
  prec_window_ = std::stod(arguments["precWindow"]);
  
  output_multiple_mass_ = (arguments["outMultipleMass"] == "true");
  
  do_final_filtering_ = (arguments["doFinalFiltering"] == "true");
  
  output_match_env_ = (arguments["outputMatchEnv"] == "true");
}

std::string DeconvPara::getParameterStr(const std::string & prefix) {
  std::stringstream output;
  output << prefix << "TopFD " << version_number << std::endl;
  // TIME_STAMP_STR is replaced later
  output << prefix << "Timestamp: " << time_util::TIME_STAMP_STR << std::endl;
  output << prefix << "********************** Parameters **********************" << std::endl;
  // output << prefix << std::setw(40) << std::left << "Input file: " << data_file_name_ << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Data type: " << "centroid" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Maximum charge: " << max_charge_ << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Maximum monoisotopic mass: " << max_mass_ << " Dalton" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Error tolerance: " << tolerance_ << " m/z" << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "MS1 signal/noise ratio: " << ms_one_sn_ratio_ << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "MS/MS signal/noise ratio: " << ms_two_sn_ratio_ << std::endl;
  output << prefix << std::setw(40) << std::left 
      << "Precursor window size: " << prec_window_ << " m/z" << std::endl;
  //output << prefix << std::setw(40) << std::left 
  //    << "Do final filtering: " << para_ptr->do_final_filtering_ << std::endl;
  output << prefix << "********************** Parameters **********************" << std::endl;
  return output.str();
}

}  // namespace toppic
