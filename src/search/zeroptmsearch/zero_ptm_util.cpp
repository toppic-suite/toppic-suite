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

#include "seq/proteoform_factory.hpp"
#include "search/zeroptmsearch/zero_ptm_util.hpp"

namespace toppic {

namespace zero_ptm_util {

PrsmPtr getPrsmPtr(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                   ZpFastMatchPtr fast_match_ptr, SpParaPtr sp_para_ptr) { 
  ProteoformPtr proteoform_ptr 
    = proteoform_factory::geneSubProteoform(fast_match_ptr->getProteoformPtr(), 
                                            fast_match_ptr->getProteoformPtr()->getFastaSeqPtr(),
                                            fast_match_ptr->getBegin(), 
                                            fast_match_ptr->getEnd());
  return getPrsmPtr(deconv_ms_ptr_vec, proteoform_ptr, sp_para_ptr);
}

PrsmPtr getPrsmPtr(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                   ProteoformPtr proteoform_ptr, SpParaPtr sp_para_ptr) {
  double refine_prec_mass = proteoform_ptr->getResSeqPtr()->getSeqMass();
  return std::make_shared<Prsm>(proteoform_ptr, deconv_ms_ptr_vec,
                                refine_prec_mass, sp_para_ptr);

}

}

}  // namespace toppic
