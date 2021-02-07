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

namespace toppic {

std::string TopfdPara::getParaStr(const std::string &prefix) {
  std::stringstream output;
  int gap = 45;
  output << prefix << "TopFD " << Version::getVersion() << std::endl;
  output << prefix << "Timestamp: " << time_util::getTimeStr() << std::endl;
  output << prefix << "###################### Parameters ######################" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Spectral data type:       " << "\t"  << "Centroid" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Maximum charge:           "  << "\t" << max_charge_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Maximum monoisotopic mass:" << 
      "\t" << max_mass_ << " Dalton" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Peak error tolerance: " << "\t" << mz_error_ << " m/z" << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "MS1 signal/noise ratio: " << "\t" << ms_one_sn_ratio_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "MS/MS signal/noise ratio: " << "\t" << ms_two_sn_ratio_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Thread number: " << "\t" << thread_number_ << std::endl;
  output << prefix << std::setw(gap) << std::left 
      << "Precursor window size: " << "\t" << prec_window_ << " m/z" << std::endl;
  if (use_env_cnn_) {
    output << prefix << std::setw(gap) << std::left 
        << "Use Env CNN model: " << "\t" << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
        << "Use Env CNN model: " << "\t" << "No" << std::endl;
  }
  if (missing_level_one_) {
    output << prefix << std::setw(gap) << std::left 
        << "Miss MS1 spectra: " << "\t" << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
        << "Miss MS1 spectra: " << "\t" << "No" << std::endl;
  }
  if (gene_html_folder_) {
    output << prefix << std::setw(gap) << std::left 
        << "Generate Html files: " << "\t" << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
        << "Generate Html files: " << "\t" << "No" << std::endl;
  }
  if (do_final_filtering_) {
    output << prefix << std::setw(gap) << std::left 
      << "Do final filtering:   " << "\t" << "Yes" << std::endl;
  }
  else {
    output << prefix << std::setw(gap) << std::left 
      << "Do final filtering:   " << "\t" << "No" << std::endl;
  }
  //output << prefix << std::setw(gap) << std::left 
     // << "Do final filtering: " << para_ptr->do_final_filtering_ << std::endl;
  output << prefix << "###################### Parameters ######################" << std::endl;
  return output.str();
}

}  // namespace toppic
