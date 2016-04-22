#include <algorithm>
#include <iostream>
#include "base/proteoform_util.hpp"
#include "zeroptmfilter/filter_protein.hpp"
#include "zeroptmfilter/mass_zero_ptm_filter.hpp"

namespace prot {

MassZeroPtmFilter::MassZeroPtmFilter(const ProteoformPtrVec &proteo_ptrs,
                                     ZeroPtmFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  std::vector<std::vector<int>> shift_2d 
      = ProteoformUtil::getNTermShift2D(proteo_ptrs, 
                                        mng_ptr->prsm_para_ptr_->getProtModPtrVec(),
                                        mng_ptr->filter_scale_);
  index_ptr_ = MassMatchPtr(new MassMatch(proteo_ptrs, shift_2d,
                                          mng_ptr->filter_scale_,
                                          mng_ptr->max_proteoform_mass_));
}

/*
void ZeroPtmFilter::computeBestMatch(const ExtendMsPtrVec &ms_ptr_vec){
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  bool pref = true;
  std::vector<std::pair<int,int>> pref_mass_errors 
      = ExtendMs::getExtendIntMassErrorList(ms_ptr_vec, pref, mng_ptr_->filter_scale_);
  pref = false;
  std::vector<std::pair<int,int>> suff_mass_errors 
      = ExtendMs::getExtendIntMassErrorList(ms_ptr_vec, pref, mng_ptr_->filter_scale_);
  std::pair<int,int> prec_minus_water_mass_error 
      = ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWaterError(tole_ptr->getPpo(), mng_ptr_->filter_scale_);
  //LOG_DEBUG("start convolution");
  index_ptr_->compZeroPtmConvolution(pref_mass_errors, suff_mass_errors, prec_minus_water_mass_error, 
                                     mng_ptr_->comp_num_, mng_ptr_->pref_suff_num_, mng_ptr_->inte_num_);

  FilterProteinPtrVec comp_prots = index_ptr_->getTopCompProts();
  comp_match_ptrs_.clear();
  int group_spec_num = ms_ptr_vec.size();
  for (size_t i = 0; i < comp_prots.size(); i++) {
    int id = comp_prots[i]->getProteinId();
    int score = comp_prots[i]->getScore();
    comp_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[id], score)));
  }

  FilterProteinPtrVec pref_prots = index_ptr_->getTopPrefProts();
  pref_match_ptrs_.clear();
  for (size_t i = 0; i < pref_prots.size(); i++) {
    int id = pref_prots[i]->getProteinId();
    int score = pref_prots[i]->getScore();
    pref_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[id], score)));
  }

  FilterProteinPtrVec suff_prots = index_ptr_->getTopSuffProts();
  suff_match_ptrs_.clear();
  for (size_t i = 0; i < suff_prots.size(); i++) {
    int id = suff_prots[i]->getProteinId();
    int score = suff_prots[i]->getScore();
    suff_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[id], score)));
  }

  FilterProteinPtrVec internal_prots = index_ptr_->getTopInternalProts();
  internal_match_ptrs_.clear();
  for (size_t i = 0; i < internal_prots.size(); i++) {
    int id = internal_prots[i]->getProteinId();
    int score = internal_prots[i]->getScore();
    internal_match_ptrs_.push_back( 
        SimplePrsmPtr(new SimplePrsm(ms_ptr_vec[0]->getMsHeaderPtr(), group_spec_num,
                                     proteo_ptrs_[id], score)));
  }

}
*/

} /* namespace prot */
