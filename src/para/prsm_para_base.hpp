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

#ifndef TOPPIC_PARA_PRSM_PARA_BASE_HPP_
#define TOPPIC_PARA_PRSM_PARA_BASE_HPP_

#include <unordered_map>

#include "para/prsm_para.hpp"

namespace toppic {

class PrsmParaBase {
 public:
  static void initBase(std::string name, 
                       std::map<std::string,std::string> &arguments);

  static PrsmParaPtr getMainPrsmParaPtr() {return main_prsm_para_ptr_;}

  static PrsmParaPtr getPrsmParaPtrByName(const std::string &name);

 private:
  static PrsmParaPtrVec prsm_para_ptr_vec_;

  static std::unordered_map<std::string, PrsmParaPtr> prsm_para_name_map_;

  static PrsmParaPtr main_prsm_para_ptr_;
};

}  // namespace toppic

#endif
