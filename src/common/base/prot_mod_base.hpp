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


#ifndef TOPPIC_COMMON_BASE_PROT_MOD_BASE_HPP_
#define TOPPIC_COMMON_BASE_PROT_MOD_BASE_HPP_

#include "common/base/ptm.hpp"
#include "common/base/trunc.hpp"
#include "common/base/prot_mod.hpp"

namespace toppic {

class ProtModBase {
 public:
  static void initBase();

  static const ProtModPtrVec& getBaseProtModPtrVec() {return prot_mod_ptr_vec_;}

  static ProtModPtr getProtModPtrByName(const std::string &name);

  static ProtModPtrVec getProtModPtrByType(const std::string &type);

  static ProtModPtr getProtModPtr_NONE() {return prot_mod_ptr_NONE_;}

  static ProtModPtr getProtModPtrFromXml(XmlDOMElement * element);

  static std::string getType_NME() {return "NME";}

  static std::string getType_NME_ACETYLATION() {return "NME_ACETYLATION";}

  static std::string getType_M_ACETYLATION() {return "M_ACETYLATION";}

  static ProtModPtr getProtModPtr_M_ACETYLATION() {return prot_mod_ptr_M_ACETYLATION_;}

 private:
  static ProtModPtrVec prot_mod_ptr_vec_;

  static ProtModPtr prot_mod_ptr_NONE_;

  static ProtModPtr prot_mod_ptr_M_ACETYLATION_;

  static std::string getName_NONE() {return "NONE";}

  static std::string getName_M_ACETYLATION() {return "M_ACETYLATION";}
};

}  // namespace toppic

#endif
