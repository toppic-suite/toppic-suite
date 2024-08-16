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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/util/time_util.hpp"

#include "common/base/base_data.hpp"
#include "ms/env/env_base.hpp"
#include "ms/mzml/pw_ms_reader.hpp"
#include "ms/mzml/mzml_profile.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"
#include "topfd/ecscore/score/onnx_ecscore.hpp"
#include "topfd/ecscore/env_coll/env_coll_detect.hpp"
#include "topfd/deconv/deconv_ms1_process.hpp"
#include "topfd/deconv/deconv_ms2_process.hpp"
#include "topdia/common/topdia_para.hpp"

//#include "topdia/pseudo_spec/generate_pseudo_spectrum.hpp"

namespace toppic {

namespace topdia_process {

void processOneFileWithFaims(TopdiaParaPtr topdia_para_ptr) {
  ////////  print parameter for each file
  /*
  std::cout << topdia_para_ptr->getParaStr("", " ");
  TopfdParaPtr para_ptr = std::make_shared<TopfdPara>(topdia_para_ptr);
  if (!para_ptr->isMissingLevelOne()) {
    std::cout << "MS1 deconvolution started." << std::endl;
    DeconvMs1ProcessPtr ms1_proc_ptr =
        std::make_shared<DeconvMs1Process>(para_ptr);
    ms1_proc_ptr->process();
    ms1_proc_ptr = nullptr;
    std::cout << "MS1 deconvolution finished." << std::endl;

    std::cout << "MS1 feature detection started." << std::endl;
    env_coll_detect::process_ms1(topdia_para_ptr);
    std::cout << "MS1 feature detection finished." << std::endl;
  }

  std::cout << "MS/MS deconvolution started." << std::endl;
  para_ptr->setMissingLevelOne(true);
  DeconvMs2ProcessPtr ms2_proc_ptr =
      std::make_shared<DeconvMs2Process>(para_ptr);
  ms2_proc_ptr->process();
  ms2_proc_ptr = nullptr;
  std::cout << "MS/MS deconvolution finished." << std::endl;

  std::cout << "MS/MS feature detection started." << std::endl;
  env_coll_detect::process_ms2(topdia_para_ptr);
  std::cout << "MS/MS feature detection finished." << std::endl;

  std::cout << "Pseudo spectrum generation started." << std::endl;
  GeneratePseudoSpectrumPtr pseudo_genrator_ptr = std::make_shared<GeneratePseudoSpectrum>(topdia_para_ptr);
  pseudo_genrator_ptr->process(topdia_para_ptr);
  pseudo_genrator_ptr = nullptr;
  std::cout << std::endl << "Pseudo spectrum generation finished." << std::endl;
  */
}

void processOneFile(TopdiaParaPtr para_ptr,
                    std::string &spec_file_name) {
  try {
    // Get mzml file profile
    PwMsReaderPtr reader_ptr = std::make_shared<PwMsReader>(spec_file_name);
    MzmlProfilePtr profile_ptr = reader_ptr->readProfile();
    int frac_id = 0;

    // check if it is faims or not
    if (profile_ptr->isFaims()) {
      bool is_faims = true;
      std::map<double, std::pair<int,int>> volt_map = profile_ptr->getVoltageMap();
      std::cout << spec_file_name << " is FAIMS data with " << volt_map.size() << " voltage levels." << std::endl;
      for (auto v : volt_map) {
        double volt = v.first;
        para_ptr->setMzmlFileNameAndFaims(spec_file_name, is_faims, volt);
        std::cout << "Processing " << spec_file_name << " with voltage " << volt << " started." << std::endl;
        para_ptr->setFracId(frac_id);
        para_ptr->setMs1ScanNumber(v.second.first);
        para_ptr->setMs2ScanNumber(v.second.second);
        processOneFileWithFaims(para_ptr);
        frac_id++;
        std::cout << "Processing " << spec_file_name << " with voltage " << volt << " finished." << std::endl;
      }
    }
    else {
      bool is_faims = false;
      double volt = -1;
      para_ptr->setMzmlFileNameAndFaims(spec_file_name, is_faims, volt);
      para_ptr->setFracId(frac_id);
      para_ptr->setMs1ScanNumber(profile_ptr->getMs1Cnt());
      para_ptr->setMs2ScanNumber(profile_ptr->getMs2Cnt());
      processOneFileWithFaims(para_ptr);
    }
  } catch (const char* e) {
    LOG_ERROR("[Exception] " << e);
    exit(EXIT_FAILURE);
  }
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

int process(TopdiaParaPtr para_ptr,  std::vector<std::string> spec_file_list) {
  // init data, envelope base, envcnn model, and ecscore model
  base_data::init();
  EnvBase::initBase(para_ptr->getResourceDir());
  onnx_env_cnn::initModel(para_ptr->getResourceDir(), para_ptr->getThreadNum());
  onnx_ecscore::initModel(para_ptr->getResourceDir(), para_ptr->getThreadNum());

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
  std::cout << "TopDIA finished." << std::endl << std::flush;
  return 0;
}

} // namespace topdia_process

}  // namespace toppic
