//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#ifndef TOPPIC_COMMON_BASE_RESIDUE_BASE_HPP_
#define TOPPIC_COMMON_BASE_RESIDUE_BASE_HPP_

#include "common/base/residue.hpp"

namespace toppic {

class ResidueBase {
 public:
  static void initBase();

  static const ResiduePtrVec& getBaseResiduePtrVec() {return residue_ptr_vec_;}

  static ResiduePtr getEmptyResiduePtr() {return empty_residue_ptr_;}

  static ResiduePtr getResiduePtrFromXml(XmlDOMElement * element);

  static ResiduePtr getBaseResiduePtr(ResiduePtr residue_ptr);

  static ResiduePtr getBaseResiduePtr(AminoAcidPtr acid_ptr, PtmPtr ptm_ptr);

  static ResiduePtr getBaseResiduePtr(AminoAcidPtr acid_ptr);

  static ResiduePtrVec getBaseNonePtmResiduePtrVec();

 private:
  static ResiduePtrVec residue_ptr_vec_;
  static ResiduePtr empty_residue_ptr_;
};

}  // namespace toppic

#endif
