//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_SPEC_PRM_MS_FACTORY_HPP_
#define PROT_SPEC_PRM_MS_FACTORY_HPP_

#include "spec/sp_para.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/prm_ms.hpp"

namespace prot {

class PrmMsFactory {
 public:
  static PrmMsPtrVec geneMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                     SpParaPtr sp_para_ptr,
                                     double prec_mono_mass);

  static PrmMsPtrVec geneSuffixMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                           SpParaPtr sp_para_ptr,
                                           double prec_mono_mass);

  static PrmMsPtrVec geneMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                     SpParaPtr sp_para_ptr,
                                     double prec_mono_mass);

  static PrmMsPtrVec geneShiftMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                          SpParaPtr sp_para_ptr, 
                                          double prec_mono_mass, double shift);
};

} /* namespace prot */

#endif /* PRM_PEAK_HPP_ */
