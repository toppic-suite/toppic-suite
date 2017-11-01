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

#include <algorithm>
#include <iostream>
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "zeroptmfilter/filter_protein.hpp"
#include "diagfilter/diag_filter.hpp"

namespace prot {

DiagFilter::DiagFilter(const ProteoformPtrVec &proteo_ptrs,
                       DiagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  bool use_rev = false;
  index_ptr_ = CompShiftPtr(new CompShift(proteo_ptrs, 
                                          mng_ptr->filter_scale_,
                                          mng_ptr->max_proteoform_mass_,
                                          mng_ptr->prsm_para_ptr_->getProtModPtrVec(),
                                          use_rev));
}

SimplePrsmPtrVec DiagFilter::getBestMatch(const PrmMsPtrVec &ms_ptr_vec){
  SimplePrsmPtrVec match_ptrs = compute(ms_ptr_vec);
  SimplePrsmPtrVec unique_match_ptrs = SimplePrsmUtil::getUniqueMatches(match_ptrs);
  std::sort(unique_match_ptrs.begin(), unique_match_ptrs.end(), SimplePrsm::cmpScoreDec);
  size_t num = mng_ptr_->filter_result_num_;
  if(num > unique_match_ptrs.size()){
    num = unique_match_ptrs.size();
  }
  SimplePrsmPtrVec result_ptrs;
  for(size_t i=0;i<num;i++){
    SimplePrsmPtr match_ptr = unique_match_ptrs[i];
    if(match_ptr->getScore() > 0.0){
      result_ptrs.push_back(match_ptr);
    }
    else{
      break;
    }
  }
  return result_ptrs;
}


inline SimplePrsmPtrVec DiagFilter::compute(const PrmMsPtrVec &ms_ptr_vec){

  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int,int>> mass_errors 
      = PrmMs::getIntMassErrorList(ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_,true,false);
  //LOG_DEBUG("mass error size " << mass_errors.size() << " filter result number " << mng_ptr_->filter_result_num_);
  SimplePrsmPtrVec match_ptrs;
  for(size_t i=0;i<mass_errors.size();i++){
    index_ptr_->compDiagConvolution(mass_errors, i, mng_ptr_->filter_result_num_);
    FilterProteinPtrVec results = index_ptr_->getTopDiagProts();
    //LOG_DEBUG("result size " << results.size());
    for(size_t j =0;j <results.size();j++){
      int id = results[j]->getProteinId();
      int score = results[j]->getScore();
      match_ptrs.push_back(
          SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(),
                                       ms_ptr_vec.size(),
                                       proteo_ptrs_[id], score)));
    }
  }
  return match_ptrs;
}

} /* namespace prot */
