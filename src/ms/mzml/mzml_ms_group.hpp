//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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
#include <vector>

#include "ms/mzml/mzml_ms.hpp"

namespace toppic {

class MzmlMsGroup {
 public:
  MzmlMsGroup(const MzmlMsPtr &ms1_ptr, const MzmlMsPtrVec &ms_ptr_vec_);

  MzmlMsPtr getMsOnePtr() const {return ms1_ptr_;}

  const MzmlMsPtrVec& getMsTwoPtrVec() const {return ms2_ptr_vec_;}

 private:
  MzmlMsPtr ms1_ptr_;
  MzmlMsPtrVec ms2_ptr_vec_;
};

using MzmlMsGroupPtr = std::shared_ptr<MzmlMsGroup>;

}  // namespace toppic

#endif 
