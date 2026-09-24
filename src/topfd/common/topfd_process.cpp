// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include "topfd/common/topfd_process.hpp"

#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "common/base/base_data.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/time_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/mzml/mzml_profile.hpp"
#include "ms/mzml/pw_ms_reader.hpp"
#include "topfd/deconv/deconv_ms1_process.hpp"
#include "topfd/deconv/deconv_ms2_process.hpp"
#include "topfd/deconv/flashdeconv_ms2_process.hpp"
#include "topfd/ecscore/env_coll/env_coll_detect.hpp"
#include "topfd/ecscore/score/onnx_ecscore.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

namespace topfd_process {

namespace {

// Intra-op threads of the ONNX sessions used by the deconvolution stages
// (see process() and processFraction()).
constexpr int kDeconvOnnxThreadNum = 1;

// Deconvolute one fraction (a single FAIMS voltage level, or the whole file
// when the data is not FAIMS): optional MS1 deconvolution and feature
// detection, followed by MS/MS deconvolution.
void processFraction(const TopfdParaPtr& para_ptr) {
  // print parameters for each fraction
  std::cout << para_ptr->getParaStr("", " ");

  if (!para_ptr->isMissingLevelOne()) {
    // MS1 spectra are present: deconvolute them and detect features first.
    std::cout << "MS1 deconvolution started." << std::endl;
    DeconvMs1ProcessPtr ms1_proc_ptr =
        std::make_shared<DeconvMs1Process>(para_ptr);
    ms1_proc_ptr->process();
    ms1_proc_ptr = nullptr;
    std::cout << "MS1 deconvolution finished." << std::endl;
    std::cout << "MS1 feature detection started." << std::endl;
    // The feature detection scores one envelope per EnvCNN call, 87,000
    // calls on 480_F3_02, each streaming the model's 87 MB of weights
    // through the CPU: memory-bound. Two intra-op threads that sleep
    // between the calls stream them fastest (st_1 feature stage: 26.1 s
    // with 1 thread, 18.7 s with 2, 19.4 s with 4, 22.9 s with 8, 27.6 s
    // with 16; toppic_test experiment 36). The single-thread session is
    // restored for the MS/MS deconvolution, whose workers call in parallel.
    int feat_onnx_thread_num = 2;
    bool feat_onnx_spinning = false;
    onnx_env_cnn::initModel(para_ptr->getResourceDir(), feat_onnx_thread_num,
                            feat_onnx_spinning);
    env_coll_detect::processMs1(para_ptr);
    onnx_env_cnn::initModel(para_ptr->getResourceDir(), kDeconvOnnxThreadNum);
    std::cout << "MS1 feature detection finished." << std::endl;
  }

  std::cout << "MS/MS deconvolution started." << std::endl;
  if (para_ptr->isUseFlashDeconv()) {
    // hidden option --flashdeconv: the FLASHDeconv algorithm
    FlashDeconvMs2ProcessPtr ms2_proc_ptr =
        std::make_shared<FlashDeconvMs2Process>(para_ptr, "ms2.msalign");
    ms2_proc_ptr->process();
  } else {
    DeconvMs2ProcessPtr ms2_proc_ptr =
        std::make_shared<DeconvMs2Process>(para_ptr, "ms2.msalign");
    ms2_proc_ptr->process();
  }
  std::cout << "MS/MS deconvolution finished." << std::endl;
}

void processOneFile(const TopfdParaPtr& para_ptr,
                    const std::string& spec_file_name) {
  try {
    // Get mzML file profile.
    PwMsReaderPtr reader_ptr = std::make_shared<PwMsReader>(spec_file_name);
    MzmlProfilePtr profile_ptr = reader_ptr->readProfile();
    para_ptr->setFilePrecWindow(profile_ptr->hasPrecWindow());

    if (profile_ptr->isFaims()) {
      // FAIMS data: process each voltage level as a separate fraction.
      const std::map<double, std::pair<int, int>>& volt_map =
          profile_ptr->getVoltageMap();
      std::cout << spec_file_name << " is FAIMS data with " << volt_map.size()
                << " voltage levels." << std::endl;
      int frac_id = 0;
      for (const auto& [volt, scan_counts] : volt_map) {
        para_ptr->setFracId(frac_id);
        para_ptr->setMs1ScanNumber(scan_counts.first);
        para_ptr->setMs2ScanNumber(scan_counts.second);
        para_ptr->setMzmlFileNameAndFaims(spec_file_name, true, volt);
        std::cout << "Processing " << spec_file_name << " with voltage " << volt
                  << " started." << std::endl;
        processFraction(para_ptr);
        std::cout << "Processing " << spec_file_name << " with voltage " << volt
                  << " finished." << std::endl;
        frac_id++;
      }
    } else {
      para_ptr->setFracId(0);
      para_ptr->setMs1ScanNumber(profile_ptr->getMs1Cnt());
      para_ptr->setMs2ScanNumber(profile_ptr->getMs2Cnt());
      para_ptr->setMzmlFileNameAndFaims(spec_file_name, false, -1);
      processFraction(para_ptr);
    }
  } catch (const char* e) {
    LOG_ERROR("[Exception] " << e);
    exit(EXIT_FAILURE);
  }
}

}  // namespace

int process(const TopfdParaPtr& para_ptr,
            const std::vector<std::string>& spec_file_list) {
  // init data, envelope base, envcnn model, and ecscore model
  base_data::init(para_ptr->getResourceDir());
  EnvBase::initBase(para_ptr->getResourceDir());
  // Both ONNX sessions run with one intra-op thread: the EnvCNN model is
  // called concurrently from the deconvolution's worker threads (one call
  // per spectrum), and the ECScore model once per candidate feature from the
  // single feature-detection thread. With intra-op threads = -u, ONNX
  // Runtime's pool (16 threads per session) oversubscribed the cores during
  // the deconvolution and spun between the per-feature calls during the
  // feature detection (toppic_test experiment 36).
  onnx_env_cnn::initModel(para_ptr->getResourceDir(), kDeconvOnnxThreadNum);
  onnx_ecscore::initModel(para_ptr->getResourceDir(), kDeconvOnnxThreadNum);

  for (const std::string& spec_file_name : spec_file_list) {
    if (!file_util::isValidMzmlFile(spec_file_name)) {
      std::cout << spec_file_name << " is not a valid mass spectral file!"
                << std::endl;
      continue;
    }
    std::cout << "Processing " << spec_file_name << " started." << std::endl;
    processOneFile(para_ptr, spec_file_name);
    std::cout << "Processing " << spec_file_name << " finished." << std::endl;
    std::cout << "Timestamp: " << time_util::getTimeStr() << std::endl;
  }

  base_data::release();
  std::cout << "TopFD finished." << std::endl << std::flush;
  return 0;
}

}  // namespace topfd_process

}  // namespace toppic
