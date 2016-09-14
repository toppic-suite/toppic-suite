// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <algorithm>
#include <iostream>
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "zeroptmfilter/filter_protein.hpp"
#include "zeroptmfilter/mass_match_factory.hpp"
#include "zeroptmfilter/mass_match_util.hpp"
#include "diagfilter/mass_diag_filter.hpp"

namespace prot {

MassDiagFilter::MassDiagFilter(const ProteoformPtrVec &proteo_ptrs,
                               DiagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  index_ptr_ = MassMatchFactory::getPrmDiagMassMatchPtr(proteo_ptrs,  
                                                        mng_ptr->max_proteoform_mass_, 
                                                        mng_ptr->filter_scale_);
}

SimplePrsmPtrVec MassDiagFilter::getBestMatch(const PrmMsPtrVec &ms_ptr_vec){
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

SimplePrsmPtrVec MassDiagFilter::compute(const PrmMsPtrVec &ms_ptr_vec){
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int,int>> mass_errors 
      = PrmMs::getIntMassErrorList(ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_,true,false);
  //LOG_DEBUG("mass error size " << mass_errors.size() << " filter result number " << mng_ptr_->filter_result_num_);
  SimplePrsmPtrVec match_ptrs;
  int row_num = index_ptr_->getRowNum();
  std::vector<int> proteo_row_begins = index_ptr_->getProteoRowBegins();
  std::vector<int> proteo_row_ends = index_ptr_->getProteoRowEnds();
  int threshold = 4;
  for(size_t i=0;i<mass_errors.size();i++){
    std::vector<short> scores(row_num, 0);
    index_ptr_->compScores(mass_errors, i, -mass_errors[i].first, scores);
    FilterProteinPtrVec results 
        = MassMatchUtil::findTopProteins(scores, proteo_row_begins, proteo_row_ends, threshold,
                                         mng_ptr_->filter_result_num_);
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
