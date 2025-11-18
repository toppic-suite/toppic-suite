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

#ifndef TOPPIC_COMMON_BASE_N_TERM_MOD_BASE_HPP_
#define TOPPIC_COMMON_BASE_N_TERM_MOD_BASE_HPP_

#include "common/base/mod.hpp"

namespace toppic {

class NTermModBase {
 public:
  static void initBase(const std::string &base_dir);

  static const ModPtrVec& getBaseNTermModPtrVec() {return mod_ptr_vec_;}

  static ModPtr getBaseNTermModPtr(ModPtr mod_ptr);

  static ModPtr getNTermNoneModPtr() {return none_mod_ptr_;}

  static bool isNTermNoneModPtr(ModPtr mod_ptr) {return mod_ptr == none_mod_ptr_;}

  static ModPtr getNTermModPtrFromXml(XmlDOMElement * element);

 private:
  static ModPtrVec mod_ptr_vec_;
  static ModPtr none_mod_ptr_;
};

}

#endif

