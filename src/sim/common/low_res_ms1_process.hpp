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

#ifndef TOPPIC_SIM_LOW_RES_MS1_PROCESS_HPP_
#define TOPPIC_SIM_LOW_RES_MS1_PROCESS_HPP_

#include "topfd/common/topfd_para.hpp"

namespace toppic {

class LowResMs1Process {
 public:
  LowResMs1Process(TopfdParaPtr topfd_para_ptr);

  void process();

 private:
  TopfdParaPtr topfd_para_ptr_;
};

using LowResMs1ProcessPtr = std::shared_ptr<LowResMs1Process>;

}

#endif
