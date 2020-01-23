//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#ifndef TOPPIC_MS_SPEC_RAW_MS_GROUP_HPP_
#define TOPPIC_MS_SPEC_RAW_MS_GROUP_HPP_

#include <memory>

#include "ms/spec/raw_ms.hpp"

namespace toppic {

class RawMsGroup {
 public:
  RawMsGroup(RawMsPtr ms1_ptr, RawMsPtrVec ms_ptr_vec_);

  RawMsPtr getMsOnePtr() {return ms1_ptr_;}

  RawMsPtrVec getMsTwoPtrVec() {return ms2_ptr_vec_;}

  int getTotalScanNum() {return total_scan_num_;}
  int setTotalScanNum(int num) {total_scan_num_ = num;}

 private:
  RawMsPtr ms1_ptr_;
  RawMsPtrVec ms2_ptr_vec_;

  int total_scan_num_;
};

typedef std::shared_ptr<RawMsGroup> RawMsGroupPtr;

} /* namespace toppic */

#endif 
