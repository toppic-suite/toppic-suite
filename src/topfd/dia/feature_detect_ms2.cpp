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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/util/file_util.hpp"

#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/feature/frac_feature_reader.hpp"
#include "topfd/ecscore/env_coll/env_coll_detect.hpp"
#include "topfd/dia/feature_detect_ms2.hpp"

namespace toppic {

FeatureDetectMs2::FeatureDetectMs2 (TopfdParaPtr topfd_para_ptr) {
  topfd_para_ptr_ = topfd_para_ptr;
}

IsolationWindowPtr FeatureDetectMs2::findWindow(IsolationWindowPtrVec & window_list, 
                                                double mz_bgn, double mz_end) {
  for (size_t i = 0; i < window_list.size(); i++) {
    if (window_list[i]->isMatch(mz_bgn, mz_end)) {
      return window_list[i];
    }
  }
  return nullptr;
}

IsolationWindowPtrVec FeatureDetectMs2::clusterIsolationWindows() {
  IsolationWindowPtrVec window_list;
  MzmlMsGroupReaderPtr reader_ptr = 
    std::make_shared<MzmlMsGroupReader>(topfd_para_ptr_->getMzmlFileName(), 
                                        topfd_para_ptr_->getPrecWindowWidth(),
                                        topfd_para_ptr_->getActivation(),
                                        topfd_para_ptr_->getFracId(),
                                        topfd_para_ptr_->isFaims(), 
                                        topfd_para_ptr_->getFaimsVoltage(), 
                                        topfd_para_ptr_->isMissingLevelOne());

  MzmlMsGroupPtr ms_group_ptr = reader_ptr->getNextMsGroupPtr();
  if (ms_group_ptr == nullptr) {
    LOG_ERROR("No spectrum to read in mzML file!");
    return window_list;
  }
  while (ms_group_ptr != nullptr) {
    MzmlMsPtrVec ms_ptr_vec = ms_group_ptr->getMsTwoPtrVec();
    for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
      MsHeaderPtr header_ptr = ms_ptr_vec[i]->getMsHeaderPtr();
      IsolationWindowPtr win_ptr = findWindow(window_list,
                                              header_ptr->getPrecWinBegin(), 
                                              header_ptr->getPrecWinEnd());
      if (win_ptr != nullptr) {
        win_ptr->addSpecId(header_ptr->getSpecId());
      }
      else {
        IsolationWindowPtr new_win_ptr =
          std::make_shared<IsolationWindow>(header_ptr->getPrecWinBegin(), 
                                            header_ptr->getPrecWinEnd());
        new_win_ptr->addSpecId(header_ptr->getSpecId());
        window_list.push_back(new_win_ptr);
      }
    }
    ms_group_ptr = reader_ptr->getNextMsGroupPtr();    
  }
  return window_list;
}

void FeatureDetectMs2::process() {
  std::string output_base_name = topfd_para_ptr_->getOutputBaseName();
  std::string ms1_feature_file_name = output_base_name + "_" + "ms1.frac_feature";
  FracFeatureReader feature_reader(ms1_feature_file_name);
  FracFeaturePtrVec ms1_feature_list = feature_reader.readAllFeatures();
  std::string ms2_pre_file_name = output_base_name + "_pre_ms2.msalign";
  std::string ms2_output_file_name = output_base_name + "_ms2.msalign";
  //file_util::rename(ms2_output_file_name, ms2_pre_file_name);
  MsAlignWriterPtr ms_writer_ptr = std::make_shared<MsAlignWriter>(ms2_output_file_name);

  IsolationWindowPtrVec window_list = clusterIsolationWindows();
  for (size_t i = 0; i < window_list.size(); i++) {
    IsolationWindowPtr win_ptr = window_list[i];
    LOG_ERROR("Mz begin " << win_ptr->getMzBgn() << " mz end " <<
              win_ptr->getMzEnd());
    env_coll_detect::processMs2(topfd_para_ptr_, ms_writer_ptr, 
                                ms1_feature_list, win_ptr->getMzBgn(),
                                win_ptr->getMzEnd(), win_ptr->getSpecIdSet());
    break;
  }
}

}

