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

#include "ms/spec/raw_ms_group.hpp"

namespace toppic {

RawMsGroup::RawMsGroup(RawMsPtr ms1_ptr, RawMsPtrVec ms2_ptr_vec):
    ms1_ptr_(ms1_ptr),
    ms2_ptr_vec_(ms2_ptr_vec) {}

} /* namespace toppic */

