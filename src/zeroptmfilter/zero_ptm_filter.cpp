#include <algorithm>
#include <iostream>
#include "zeroptmfilter/zero_ptm_filter.hpp"

namespace prot {

ZeroPtmFilter::ZeroPtmFilter(const ProteoformPtrVec &proteo_ptrs,
                             ZeroPtmFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  index_ptr_ = CompShiftPtr( new CompShift(proteo_ptrs, 
                                           mng_ptr->filter_scale_,
                                           mng_ptr->max_proteoform_mass_,
                                           mng_ptr->prsm_para_ptr_->getProtModPtrVec()));
}

void ZeroPtmFilter::computeBestMatch(const ExtendMsPtrVec &ms_ptr_vec){
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<std::pair<int,int>> mass_errors 
      = ExtendMs::getExtendIntMassErrorList(ms_ptr_vec, mng_ptr_->filter_scale_);
  std::pair<int,int> prec_mass_error 
      = ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWaterError(mng_ptr_->filter_scale_);
  //LOG_DEBUG("start convolution");
  index_ptr_->compZeroPtmConvolution(mass_errors, prec_mass_error, mng_ptr_->comp_num_, 
                                     mng_ptr_->pref_suff_num_, mng_ptr_->inte_num_);

  std::vector<std::pair<int,int>> comp_scores = index_ptr_->getTopCompProteoScores();
  comp_match_ptrs_.clear();
  int group_spec_num = ms_ptr_vec.size();
  for (size_t i = 0; i < comp_scores.size(); i++) {
    comp_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[comp_scores[i].first], comp_scores[i].second)));
  }

  std::vector<std::pair<int,int>> pref_scores = index_ptr_->getTopPrefProteoScores();
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_scores.size(); i++) {
    pref_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[pref_scores[i].first], pref_scores[i].second)));
  }

  std::vector<std::pair<int,int>> suff_scores = index_ptr_->getTopSuffProteoScores();
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_scores.size(); i++) {
    suff_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[suff_scores[i].first], suff_scores[i].second)));
  }

  std::vector<std::pair<int,int>> internal_scores = index_ptr_->getTopInternalProteoScores();
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_scores.size(); i++) {
    internal_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[internal_scores[i].first], internal_scores[i].second)));
  }

}

} /* namespace prot */
