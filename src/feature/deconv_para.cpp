//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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
#include <map>
#include <string>

#include "feature/deconv_para.hpp"

namespace prot {

DeconvPara::DeconvPara(std::map<std::string, std::string> &arguments) {
  data_file_name_ = arguments["spectrumFileName"];

  missing_level_one_ = (arguments["missingLevelOne"] == "true");
  max_charge_ = std::stoi(arguments["maxCharge"]);
  max_mass_ = std::stod(arguments["maxMass"]);
  tolerance_ = std::stod(arguments["mzError"]);
  ms_two_sn_ratio_ = std::stod(arguments["msTwoSnRatio"]);
  ms_one_sn_ratio_ = std::stod(arguments["msOneSnRatio"]);
  keep_unused_peaks_ = (arguments["keepUnusedPeaks"] == "true");
  prec_window_ = std::stod(arguments["precWindow"]);
  exec_dir_ = arguments["executiveDir"];
  output_multiple_mass_ = (arguments["outMultipleMass"] == "true");
  do_final_filtering_ = (arguments["doFinalFiltering"] == "true");
  output_match_env_ = (arguments["outputMatchEnv"] == "true");
  //std::cout << "Do final filtering " << do_final_filtering_;
}

}  // namespace prot
