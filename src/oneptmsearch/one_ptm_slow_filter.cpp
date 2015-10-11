#include "oneptmsearch/one_ptm_slow_filter.hpp"

namespace prot {

OnePtmSlowFilter::OnePtmSlowFilter(
    SpectrumSetPtr spectrum_set_ptr,
    SimplePrsmPtrVec simple_prsm_ptrs,
    CompShiftLowMemPtr comp_shift_ptr,
    SemiAlignTypePtr type_ptr, 
    OnePtmSearchMngPtr mng_ptr){

  OnePtmSlowMatchPtrVec slow_match_ptrs;
  // init complete_prefix_slow_match_ptrs
  for(size_t i=0;i<simple_prsm_ptrs.size();i++){
    if (type_ptr == SemiAlignTypeFactory::getCompletePtr() 
        || type_ptr == SemiAlignTypeFactory::getPrefixPtr()) {
      ProteoformPtrVec mod_proteo_ptrs = simple_prsm_ptrs[i]->getModProteoformPtrs(); 
      for (size_t j = 0; j < mod_proteo_ptrs.size(); j++) {
        OnePtmSlowMatchPtr one_ptm_slow_match_ptr(
            new OnePtmSlowMatch(mod_proteo_ptrs[j],spectrum_set_ptr, comp_shift_ptr, type_ptr, mng_ptr));
        slow_match_ptrs.push_back(one_ptm_slow_match_ptr);
      }
    }
    else {
      ProteoformPtr raw_proteo_ptr = simple_prsm_ptrs[i]->getProteoformPtr();
      OnePtmSlowMatchPtr one_ptm_slow_match_ptr(
          new OnePtmSlowMatch(raw_proteo_ptr, spectrum_set_ptr, comp_shift_ptr, type_ptr, mng_ptr));
      slow_match_ptrs.push_back(one_ptm_slow_match_ptr);
    }
  }

  for(size_t i=0; i<slow_match_ptrs.size();i++){
    PrsmPtr prsm_ptr = slow_match_ptrs[i]->getResultPrsm();
    prsm_ptrs_.push_back(prsm_ptr);
  }
}

} /* namespace prot */
