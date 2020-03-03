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


#ifndef TOPPIC_MS_SPEC_EXTEND_MS_FACTORY_HPP_
#define TOPPIC_MS_SPEC_EXTEND_MS_FACTORY_HPP_

#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/extend_ms.hpp"
#include "para/sp_para.hpp"

namespace toppic {

namespace extend_ms_factory {

ExtendMsPtr geneMsThreePtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr,
                           double new_prec_mass);

ExtendMsPtrVec geneMsThreePtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                 SpParaPtr sp_para_ptr, double new_prec_mass);

}  // namespace extend_ms_factory

}  // namespace toppic

#endif 
