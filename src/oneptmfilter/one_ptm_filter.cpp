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
#include "zeroptmfilter/filter_protein.hpp"
#include "oneptmfilter/one_ptm_filter.hpp"

namespace prot {

OnePtmFilter::OnePtmFilter(const ProteoformPtrVec &proteo_ptrs,
                           OnePtmFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  bool use_rev = true;
  index_ptr_ = CompShiftPtr(new CompShift(proteo_ptrs, 
                                          mng_ptr->filter_scale_,
                                          mng_ptr->max_proteoform_mass_,
                                          mng_ptr->prsm_para_ptr_->getProtModPtrVec(),
                                          use_rev));
}

void OnePtmFilter::computeBestMatch(const PrmMsPtrVec &prm_ms_ptr_vec, 
                                    const PrmMsPtrVec &srm_ms_ptr_vec){
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int,int>> pref_mass_errors 
      = PrmMs::getIntMassErrorList(prm_ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_, true, false);
  std::vector<std::pair<int,int>> suff_mass_errors 
      = PrmMs::getIntMassErrorList(srm_ms_ptr_vec, tole_ptr, mng_ptr_->filter_scale_, false, true);
  //LOG_DEBUG("start convolution");
  index_ptr_->compOnePtmConvolution(pref_mass_errors, suff_mass_errors, 
                                    mng_ptr_->comp_num_, mng_ptr_->pref_suff_num_, 
                                    mng_ptr_->inte_num_, mng_ptr_->shift_num_);
  FilterProteinPtrVec comp_prots = index_ptr_->getTopCompProts();
  comp_match_ptrs_.clear();
  int group_spec_num = prm_ms_ptr_vec.size();
  for (size_t i = 0; i < comp_prots.size(); i++) {
    int id = comp_prots[i]->getProteinId();
    int score = comp_prots[i]->getScore();
    comp_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(prm_ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[id], score)));
  }

  FilterProteinPtrVec pref_prots = index_ptr_->getTopPrefProts();
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_prots.size(); i++) {
    int id = pref_prots[i]->getProteinId();
    int score = pref_prots[i]->getScore();
    SimplePrsmPtr prsm_ptr(new SimplePrsm(prm_ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                          proteo_ptrs_[id], score));
    prsm_ptr->setCTruncShifts(pref_prots[i]->getCTermShifts());
    pref_match_ptrs_.push_back(prsm_ptr); 
  }

  FilterProteinPtrVec suff_prots = index_ptr_->getTopSuffProts();
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_prots.size(); i++) {
    int id = suff_prots[i]->getProteinId();
    int score = suff_prots[i]->getScore();
    SimplePrsmPtr prsm_ptr(new SimplePrsm(prm_ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                          proteo_ptrs_[id], score));
    prsm_ptr->setNTruncShifts(suff_prots[i]->getNTermShifts());
    suff_match_ptrs_.push_back(prsm_ptr); 
  }

  FilterProteinPtrVec internal_prots = index_ptr_->getTopInternalProts();
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_prots.size(); i++) {
    int id = internal_prots[i]->getProteinId();
    int score = internal_prots[i]->getScore();
    SimplePrsmPtr prsm_ptr(new SimplePrsm(prm_ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                          proteo_ptrs_[id], score));
    prsm_ptr->setNTruncShifts(internal_prots[i]->getNTermShifts());
    prsm_ptr->setCTruncShifts(internal_prots[i]->getCTermShifts());
    internal_match_ptrs_.push_back(prsm_ptr); 
  }
}

} /* namespace prot */
