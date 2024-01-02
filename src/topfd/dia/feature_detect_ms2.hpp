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

#ifndef TOPPIC_TOPFD_DIA_FEATURE_DETECT_MS2_HPP_
#define TOPPIC_TOPFD_DIA_FEATURE_DETECT_MS2_HPP_

#include "topfd/common/topfd_para.hpp"

#include "topfd/dia/isolation_window.hpp"

namespace toppic {

class FeatureDetectMs2 {
 public:
  FeatureDetectMs2(TopfdParaPtr topfd_para_ptr);

  void process();

 private:
  IsolationWindowPtrVec clusterIsolationWindows(); 

  IsolationWindowPtr findWindow(IsolationWindowPtrVec & window_list, 
                                double mz_bgn, double mz_end);

  void processSingleWindow(IsolationWindowPtr win_ptr);

  TopfdParaPtr topfd_para_ptr_;
  
};

  
typedef std::shared_ptr<FeatureDetectMs2> FeatureDetectMs2Ptr;

}

#endif
