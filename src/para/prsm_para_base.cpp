//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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
#include "para/prsm_para_base.hpp"

namespace toppic {

PrsmParaPtrVec PrsmParaBase::prsm_para_ptr_vec_;

PrsmParaPtr PrsmParaBase::main_prsm_para_ptr_;

std::unordered_map<std::string, PrsmParaPtr> PrsmParaBase::prsm_para_name_map_;

void PrsmParaBase::initBase(std::string name, 
                            std::map<std::string,std::string> &arguments) {
  PrsmParaPtr ptr_ = std::make_shared<PrsmPara>(name, arguments);

  main_prsm_para_ptr_ = ptr_;

  prsm_para_ptr_vec_.push_back(ptr_);

  prsm_para_name_map_[ptr_->getName()] = ptr_;
}

PrsmParaPtr PrsmParaBase::getPrsmParaPtrByName(const std::string &name) {
  PrsmParaPtr prsm_para_ptr = prsm_para_name_map_[name];
  if (prsm_para_ptr == nullptr) {
    LOG_ERROR("Sp para " << name << " cannot be found!");
  }
  return prsm_para_ptr;
}

}  // namespace toppic
