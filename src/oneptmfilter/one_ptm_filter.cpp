#include <algorithm>
#include <iostream>
#include "oneptmfilter/one_ptm_filter.hpp"

namespace prot {

OnePtmFilter::OnePtmFilter(const ProteoformPtrVec &proteo_ptrs,
                           OnePtmFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  index_ptr_ = OnePtmCompShiftPtr(new OnePtmCompShift(proteo_ptrs, mng_ptr));
}

/*
void OnePtmFilter::compBestMatch(PrmMsPtr ms_ptr){
  SimplePrsmPtrVec match_ptrs = compute(ms_ptr);
  SimplePrsmPtrVec unique_match_ptrs = getUniqueMatches(match_ptrs);
  std::sort(unique_match_ptrs.begin(), unique_match_ptrs.end(),simplePrsmDown);
  size_t num = mng_ptr_->one_ptm_filter_result_num_;
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
*/

void OnePtmFilter::computeBestMatch(PrmMsPtr ms_ptr){
  std::pair<std::vector<int>, std::vector<int>> mass_errors 
      = getIntMassErrorList(ms_ptr, mng_ptr_->ptm_fast_filter_scale_, true, false);
  //LOG_DEBUG("start convolution");
  index_ptr_->compConvolution(mass_errors.first, mass_errors.second, 
                              mng_ptr_->one_ptm_filter_result_num_);
  std::vector<std::pair<int,int>> comp_scores = index_ptr_->getTopCompProteoScores();
  comp_match_ptrs_.clear();
  for (size_t i = 0; i < comp_scores.size(); i++) {
    comp_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr->getHeaderPtr(),
                                     proteo_ptrs_[comp_scores[i].first], comp_scores[i].second)));
  }

  std::vector<std::pair<int,int>> pref_scores = index_ptr_->getTopPrefProteoScores();
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_scores.size(); i++) {
    pref_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr->getHeaderPtr(),
                                     proteo_ptrs_[pref_scores[i].first], pref_scores[i].second)));
  }

  std::vector<std::pair<int,int>> suff_scores = index_ptr_->getTopSuffProteoScores();
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_scores.size(); i++) {
    suff_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr->getHeaderPtr(),
                                     proteo_ptrs_[suff_scores[i].first], suff_scores[i].second)));
  }

  std::vector<std::pair<int,int>> internal_scores = index_ptr_->getTopInternalProteoScores();
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_scores.size(); i++) {
    internal_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr->getHeaderPtr(),
                                     proteo_ptrs_[internal_scores[i].first], internal_scores[i].second)));
  }

}


} /* namespace prot */
