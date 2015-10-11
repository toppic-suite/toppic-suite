#include "ptmsearch/ptm_slow_filter.hpp"

namespace prot {

PtmSlowFilter::PtmSlowFilter(
    SpectrumSetPtr spectrum_set_ptr,
    SimplePrsmPtrVec simple_prsm_ptrs,
    CompShiftLowMemPtr comp_shift_ptr,
    PtmMngPtr mng_ptr){

  // init complete_prefix_slow_match_ptrs
  for(size_t i=0;i<simple_prsm_ptrs.size();i++){
    ProteoformPtrVec raw_proteo_ptrs;
    raw_proteo_ptrs.push_back(simple_prsm_ptrs[i]->getProteoformPtr());
    ProteoformPtrVec mod_proteo_ptrs = simple_prsm_ptrs[i]->getModProteoformPtrs(); 
    for (size_t j = 0; j < mod_proteo_ptrs.size(); j++) {
      PtmSlowMatchPtr ptm_slow_match_ptr(new PtmSlowMatch(mod_proteo_ptrs[j],spectrum_set_ptr,
                                                          comp_shift_ptr,mng_ptr));
      complete_prefix_slow_match_ptrs_.push_back(ptm_slow_match_ptr);
    }
  }

  // init suffix_internal_slow_prsm_ptrs
  for(size_t i=0; i<complete_prefix_slow_match_ptrs_.size();i++){
    ProtModPtr prot_mod_ptr = complete_prefix_slow_match_ptrs_[i]->getProteoform()->getProtModPtr();
    if (prot_mod_ptr == ProtModFactory::getProtModPtr_NONE()) {
      suffix_internal_slow_match_ptrs_.push_back(complete_prefix_slow_match_ptrs_[i]);
    }
  }

  // compute complete and prefix prsms 
  for(size_t i=0; i<complete_prefix_slow_match_ptrs_.size();i++){
    PrsmPtrVec comp_ptrs;
    //LOG_DEBUG("Compute complete prsm " << i);
    complete_prefix_slow_match_ptrs_[i]->compute(SemiAlignTypeFactory::getCompletePtr(), comp_ptrs);
    complete_prsm_2d_ptrs_.push_back(comp_ptrs);
    PrsmPtrVec prefix_ptrs;
    //LOG_DEBUG("Compute prefix prsm " << i);
    complete_prefix_slow_match_ptrs_[i]->compute(SemiAlignTypeFactory::getPrefixPtr(), prefix_ptrs);
    prefix_prsm_2d_ptrs_.push_back(prefix_ptrs);
    //LOG_DEBUG("compute prefi completed");
  }
  //LOG_DEBUG("complete prefix completed");

  // compute suffix and internal prsms 
  for(size_t i=0; i< suffix_internal_slow_match_ptrs_.size();i++){
    PrsmPtrVec suffix_ptrs;
    suffix_internal_slow_match_ptrs_[i]->compute(SemiAlignTypeFactory::getSuffixPtr(), suffix_ptrs);
    suffix_prsm_2d_ptrs_.push_back(suffix_ptrs);
    PrsmPtrVec internal_ptrs;
    suffix_internal_slow_match_ptrs_[i]->compute(SemiAlignTypeFactory::getInternalPtr(), internal_ptrs);
    internal_prsm_2d_ptrs_.push_back(internal_ptrs);
  }
  //LOG_DEBUG("suffix internal completed");
}

PrsmPtrVec PtmSlowFilter::getPrsms(int shift_num, SemiAlignTypePtr type_ptr){
  PrsmPtrVec prsm_ptrs;
  if (type_ptr == SemiAlignTypeFactory::getCompletePtr()) {
    for (size_t i = 0; i < complete_prsm_2d_ptrs_.size(); i++) {
      if (complete_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(complete_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  else if (type_ptr == SemiAlignTypeFactory::getPrefixPtr()) {
    for (size_t i = 0; i < prefix_prsm_2d_ptrs_.size(); i++) {
      if (prefix_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(prefix_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  else if (type_ptr == SemiAlignTypeFactory::getSuffixPtr()) {
    for (size_t i = 0; i < suffix_prsm_2d_ptrs_.size(); i++) {
      if (suffix_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(suffix_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  else if (type_ptr == SemiAlignTypeFactory::getInternalPtr()) {
    for (size_t i = 0; i < internal_prsm_2d_ptrs_.size(); i++) {
      if (internal_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(internal_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  return prsm_ptrs;
}

} /* namespace prot */
