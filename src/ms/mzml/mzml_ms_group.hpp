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

#ifndef TOPPIC_MS_MZML_MZML_MS_GROUP_HPP_
#define TOPPIC_MS_MZML_MZML_MS_GROUP_HPP_

#include <memory>

#include "ms/mzml/mzml_ms.hpp"

namespace toppic {

class MzmlMsGroup {
 public:
  MzmlMsGroup(MzmlMsPtr ms1_ptr, MzmlMsPtrVec ms_ptr_vec_);

  MzmlMsPtr getMsOnePtr() {return ms1_ptr_;}

  MzmlMsPtrVec getMsTwoPtrVec() {return ms2_ptr_vec_;}

 private:
  MzmlMsPtr ms1_ptr_;
  MzmlMsPtrVec ms2_ptr_vec_;
};

typedef std::shared_ptr<MzmlMsGroup> MzmlMsGroupPtr;

} /* namespace toppic */

#endif 
