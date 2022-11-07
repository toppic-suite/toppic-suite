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
#include <sstream>
#include <iomanip>
#include "common/util/version.hpp"
#include "common/util/time_util.hpp"
#include "topfd/common/topfd_para.hpp"
#include "common/util/version.hpp"

namespace toppic {
std::string TopfdPara::getParaStr(const std::string &prefix, 
		                  const std::string &sep) {
  std::stringstream output;
  int gap = 25;
  output << prefix << "TopFD " << Version::getVersion() << std::endl;
  output << prefix << "Timestamp: " << time_util::getTimeStr() << std::endl;
  output << prefix << "###################### Parameters ######################" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Activation type:          " << sep  << activation_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Number of MS1 scans:      " << sep  << ms_1_scan_num_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Number of MS/MS scans:    " << sep  << ms_2_scan_num_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Spectral data type:       " << sep  << "Centroid" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Maximum charge:           "  << sep << max_charge_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Maximum monoisotopic mass:" << sep << max_mass_ << " Dalton" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Peak error tolerance:     " << sep << mz_error_ << " m/z" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "MS1 signal/noise ratio:   " << sep << ms_one_sn_ratio_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "MS/MS signal/noise ratio: " << sep << ms_two_sn_ratio_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Thread number:            " << sep << thread_num_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Precursor window size:    " << sep << prec_window_ << " m/z" << std::endl;
  if (use_env_cnn_) {
    output << prefix << std::setw(gap) << std::left 
      << "Use Env CNN model:        " << sep << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Use Env CNN model:        " << sep << "No" << std::endl;
  }
  if (missing_level_one_) {
    output << prefix << std::setw(gap) << std::left 
      << "Miss MS1 spectra:         " << sep << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Miss MS1 spectra:         " << sep << "No" << std::endl;
  }
  if (gene_html_folder_) {
    output << prefix << std::setw(gap) << std::left 
      << "Generate Html files:      " << sep << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Generate Html files:      " << sep << "No" << std::endl;
  }
  if (do_final_filtering_) {
    output << prefix << std::setw(gap) << std::left 
      << "Do final filtering:       " << sep << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Do final filtering:       " << sep << "No" << std::endl;
  }
  output << prefix << std::setw(gap) << std::left 
      << "Version:                  " << sep << Version::getVersion() << std::endl;   
  output << prefix << "###################### Parameters ######################" << std::endl;
  return output.str();
}

}  // namespace toppic
