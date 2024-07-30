//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
//
/*
#include <iomanip>

#include "common/util/version.hpp"
#include "common/base/mass_constant.hpp"
#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/feature/feature_merge.hpp"
#include "topfd/msreader/raw_ms_writer.hpp"
#include "topfd/deconv/deconv_process.hpp"
#include "topfd/feature_detect/feature_detect.hpp"
*/

#include <fstream>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/util/time_util.hpp"
#include "common/base/base_data.hpp"

#include "ms/env/env_base.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp" 
#include "topfd/deconv/deconv_single_sp.hpp"

namespace toppic {

namespace topfd_single_process {

PeakPtrVec readPeakFile(std::string file_name) {
  PeakPtrVec peak_list;
  std::ifstream input;
  input.open(file_name.c_str(), std::ios::in);
  std::string line;
  while (std::getline(input, line)) {
    str_util::trim(line);
    if (line.length() == 0) {
      continue;
    } 
    std::vector<std::string> strs = str_util::split(line, " ");
    double mz = std::stod(strs[0]);
    double inte = std::stod(strs[1]);
    //std::cout << "mz " << mz << " intensity " << inte << std::endl;
    PeakPtr peak_ptr = std::make_shared<Peak>(mz, inte);
    peak_list.push_back(peak_ptr);
  }
  input.close();
  //std::cout << "peak list finished" << std::endl;
  return peak_list;
}

int processOneFile(TopfdParaPtr para_ptr, 
    const std::string &spec_file_name) {
  try {
    int ms_level = 2;
    double max_mass = para_ptr->getMaxMass();
    double max_charge = para_ptr->getMaxCharge();
    PeakPtrVec peak_list = readPeakFile(spec_file_name); 

    MatchEnvPtrVec result_envs; 
    if (peak_list.size() > 0) {
      DeconvSingleSpPtr deconv_ptr = std::make_shared<DeconvSingleSp>(para_ptr, peak_list, 
          ms_level, max_mass, max_charge);
      result_envs = deconv_ptr->deconv(); 
    }

    // header
    MsHeaderPtr header_ptr = std::make_shared<MsHeader>();
    header_ptr->setSpecId(0);
    header_ptr->setSingleScan(1);

    DeconvMsPtr ms_ptr = match_env_util::getDeconvMsPtr(header_ptr, result_envs);

    std::string output_base_name = para_ptr->getOutputBaseName();
    std::string ms2_msalign_name = output_base_name + "_ms2.msalign";
    MsAlignWriterPtr ms2_writer_ptr = std::make_shared<MsAlignWriter>(ms2_msalign_name);
    ms2_writer_ptr->writeMs(ms_ptr);
    ms2_writer_ptr = nullptr;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

bool isValidFile(std::string &file_name) {
  if (str_util::endsWith(file_name, "mzML")
      || str_util::endsWith(file_name, "mzXML")
      || str_util::endsWith(file_name, "mzml")
      || str_util::endsWith(file_name, "mzxml")) {
    return true;
  }
  else {
    return false;
  }
}


int process(TopfdParaPtr para_ptr,  std::vector<std::string> spec_file_list) {
  // init data, envelope base, envcnn model, and ecscore model
  base_data::init();
  EnvBase::initBase(para_ptr->getResourceDir());
  onnx_env_cnn::initModel(para_ptr->getResourceDir(), para_ptr->getThreadNum());

  for (size_t k = 0; k < spec_file_list.size(); k++) {
    if (isValidFile(spec_file_list[k])) {
      std::cout << "Processing " << spec_file_list[k] << " started." << std::endl;
      processOneFile(para_ptr, spec_file_list[k]); 
      std::cout << "Timestamp: " << time_util::getTimeStr() << std::endl;
      std::cout << "Processing " << spec_file_list[k] << " finished." << std::endl;
    }
    else {
      std::cout << spec_file_list[k] << " is not a valid mass spectral file!" << std::endl; 
    }
  }

  base_data::release();
  std::cout << "TopFD single finished." << std::endl << std::flush;
  return 0;
}



} // namespace topfd_process 

}  // namespace toppic

