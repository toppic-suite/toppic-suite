//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include <ctime>
#include <cstdlib> 
#include <iostream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/spec/msalign_reader_util.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/para/ecscore_para.hpp"
#include "sim/common/ms_feat_util.hpp"
/*
#include "ms/feature/spec_feature_writer.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/ecscore/env/seed_env.hpp"
#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"
#include "topfd/ecscore/score/ecscore.hpp"
#include "topfd/ecscore/score/ecscore_writer.hpp"
#include "topfd/ecscore/env_coll/env_coll_assign.hpp"
*/

namespace toppic {

namespace low_res_feat {

void processMs1(TopfdParaPtr topfd_para_ptr) {
  if (topfd_para_ptr->isMissingLevelOne()) {
    return;
  }
  
  // read deconvoluted MS1 peaks
  std::string output_base_name = topfd_para_ptr->getOutputBaseName();

  // read ms1 raw peaks and ms2_headers
  PeakPtrVec2D ms1_mzml_peaks;
  MsHeaderPtr2D ms2_header_ptr_2d;
  MzmlMsGroupReaderPtr mzml_reader_ptr = 
    std::make_shared<MzmlMsGroupReader>(topfd_para_ptr->getMzmlFileName(), 
                                        topfd_para_ptr->getPrecWindowWidth(),
                                        topfd_para_ptr->getActivation(),
                                        topfd_para_ptr->getFracId(),
                                        topfd_para_ptr->isFaims(), 
                                        topfd_para_ptr->getFaimsVoltage(), 
                                        topfd_para_ptr->isMissingLevelOne());
  mzml_reader_ptr->getMs1Map(ms1_mzml_peaks, ms2_header_ptr_2d); 
  MsHeaderPtrVec ms1_header_ptr_vec = mzml_reader_ptr->getMs1HeaderPtrVec();
  mzml_reader_ptr = nullptr;

  double sn_ratio = topfd_para_ptr->getMsOneSnRatio();
  bool single_scan_noise = topfd_para_ptr->isUseSingleScanNoiseLevel();

  EcscoreParaPtr score_para_ptr = std::make_shared<EcscorePara>(
      topfd_para_ptr->getFracId(), topfd_para_ptr->getMzmlFileName(),
      topfd_para_ptr->getMaxCharge(), topfd_para_ptr->getMs1MinScanNum());

  MsMapPtr matrix_ptr = std::make_shared<MsMap>(ms1_mzml_peaks, ms1_header_ptr_vec, 
                                                score_para_ptr->bin_size_,
                                                sn_ratio, single_scan_noise);

  if (score_para_ptr->min_scan_num_ >= 2) {
    matrix_ptr->removeNonNeighbors(score_para_ptr->neighbor_mz_tole_);
  }

  ms_feat_util::writepng(matrix_ptr, output_base_name + "_ms1.png");

  /*
  for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
    int count = seed_env_idx + 1;
    if (count % 100 == 0 || count == seed_num) {
      double perc = static_cast<int>(count * 100 / seed_num);
      std::cout << "\r" << "Processing feature " << count << " ...       " << perc << "\% finished." << std::flush;
    }
    env_coll_ptr->setEcscore(ecscore_ptr->getScore());
    env_coll_ptr->removePeakData(matrix_ptr);
    env_coll_list.push_back(env_coll_ptr);
  }
  */
  std::cout << "\n";
}

}  // namespace

}  // namespace toppic
