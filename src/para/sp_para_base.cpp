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
#include "common/base/mass_constant.hpp"
#include "common/base/activation_base.hpp"
#include "para/sp_para_base.hpp"

namespace toppic {

SpParaPtrVec SpParaBase::sp_para_ptr_vec_;

SpParaPtr SpParaBase::main_sp_para_ptr_;

std::unordered_map<std::string, SpParaPtr> SpParaBase::sp_para_name_map_;

void SpParaBase::initBase(std::string name, std::string activation_name) {
  double IM = mass_constant::getIsotopeMass();
  // the set of offsets used to expand the monoisotopic mass list 
  std::vector<double> ext_offsets {{0, -IM, IM}};
  double extend_min_mass = 5000;

  int min_peak_num = 10;
  double min_mass = 50.0;

  ActivationPtr activation_ptr; 
  if (activation_name != "FILE") {
    activation_ptr = ActivationBase::getActivationPtrByName(activation_name);
  }
  SpParaPtr ptr_ = std::make_shared<SpPara>(name, min_peak_num, min_mass, 
                                            extend_min_mass, ext_offsets, 
                                            activation_ptr);
  main_sp_para_ptr_ = ptr_;

  sp_para_ptr_vec_.push_back(ptr_);

  sp_para_name_map_[ptr_->getName()] = ptr_;
}


SpParaPtr SpParaBase::getSpParaPtrByName(const std::string &name) {
  SpParaPtr sp_para_ptr = sp_para_name_map_[name];
  if (sp_para_ptr == nullptr) {
    LOG_ERROR("Sp para " << name << " cannot be found!");
  }
  return sp_para_ptr;
}

}  // namespace toppic
