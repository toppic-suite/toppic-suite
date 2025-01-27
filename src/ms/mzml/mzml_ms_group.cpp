//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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
#include "ms/mzml/mzml_ms_group.hpp"

namespace toppic {

MzmlMsGroup::MzmlMsGroup(MzmlMsPtr ms1_ptr, 
                         MzmlMsPtrVec ms2_ptr_vec):
    ms1_ptr_(ms1_ptr),
    ms2_ptr_vec_(ms2_ptr_vec) {
      //init ms one id
      for (size_t i = 0; i < ms2_ptr_vec.size(); i++) {
        if (ms1_ptr_ == nullptr) {
          ms2_ptr_vec[i]->getMsHeaderPtr()->setMsOneId(-1); 
        }
        else {
          ms2_ptr_vec[i]->getMsHeaderPtr()->setMsOneId(ms1_ptr->getMsHeaderPtr()->getSpecId());
        }
      }
    }
} /* namespace toppic */

