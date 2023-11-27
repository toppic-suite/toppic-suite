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

#ifndef TOPPIC_TOPFD_DECONV_MS1_PROCESS_HPP_
#define TOPPIC_TOPFD_DECONV_MS1_PROCESS_HPP_

#include "ms/env/env_para.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/dp/dp_para.hpp"

namespace toppic {

class DeconvMs1Process {
 public:
  DeconvMs1Process(TopfdParaPtr topfd_para_ptr);

  void process();

 private:
  TopfdParaPtr topfd_para_ptr_;
  EnvParaPtr env_para_ptr_;
  DpParaPtr dp_para_ptr_;
  
  void prepareFileFolder();
};

typedef std::shared_ptr<DeconvMs1Process> DeconvMs1ProcessPtr;

}

#endif
