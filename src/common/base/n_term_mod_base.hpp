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

#include "common/base/n_term_mod.hpp"

namespace toppic {

class NTermModBase {
 public:
  static void initBase(const std::string &base_dir);

  static const NTermModPtrVec& getBaseNTermModPtrVec() {return mod_ptr_vec_;}

  static NTermModPtr getBaseNTermModPtr(NTermModPtr n_term_mod_ptr);

  static NTermModPtr getNoneModPtr() {return none_mod_ptr_;}

  static bool isNoneModPtr(NTermModPtr n_term_mod_ptr) {return n_term_mod_ptr == none_mod_ptr_;}

  static NTermModPtr getNTermModPtrFromXml(XmlDOMElement * element);

 private:
  static NTermModPtrVec mod_ptr_vec_;
  static NTermModPtr none_mod_ptr_;
};

}

#endif

