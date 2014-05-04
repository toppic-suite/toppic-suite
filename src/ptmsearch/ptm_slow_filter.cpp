/*
 * ptm_slow_filter.cpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#include "ptmsearch/ptm_slow_filter.hpp"

namespace prot {
PtmSlowFilter::PtmSlowFilter(
    SpectrumSetPtr spectrum_set,
    SimplePrSMPtrVec fast_Matches,
    CompShiftLowMemPtr comp_shift,
    PtmMngPtr mng){

  // init complete_prefix_slow_matches
  for(unsigned int i=0;i<fast_Matches.size();i++){
    ProteoformPtrVec raw_forms;
    raw_forms.push_back(fast_Matches[i]->getSeq());
    ProteoformPtrVec prot_mod_forms 
        = generateProtModProteoform(raw_forms, mng->prsm_para_ptr_->getAllowProtModPtrVec());
    for (unsigned int j = 0; j < prot_mod_forms.size(); j++) {
      complete_prefix_slow_matches_.push_back(PtmSlowMatchPtr(
            new PtmSlowMatch(prot_mod_forms[j],spectrum_set,comp_shift,mng)));
    }
  }
  // init suffix_internal_slow_matches
  for(unsigned int i=0; i<complete_prefix_slow_matches_.size();i++){
    ProtModPtr prot_mod = complete_prefix_slow_matches_[i]->getSeq()->getProtModPtr();
    if (prot_mod == ProtModFactory::getProtModPtr_NONE()) {
      suffix_internal_slow_matches_.push_back(complete_prefix_slow_matches_[i]);
    }
  }

  // compute complete and prefix prsms 
  for(unsigned int i=0; i<complete_prefix_slow_matches_.size();i++){
    PrSMPtrVec comps;
    complete_prefix_slow_matches_[i]->compute(SemiAlignTypeFactory::getCompletePtr(), comps);
    complete_prsms_.push_back(comps);
    PrSMPtrVec prefixs;
    complete_prefix_slow_matches_[i]->compute(SemiAlignTypeFactory::getPrefixPtr(), prefixs);
    prefix_prsms_.push_back(prefixs);
  }

  // compute suffix and internal prsms 
  for(unsigned int i=0; i< suffix_internal_slow_matches_.size();i++){
    PrSMPtrVec suffixs;
    suffix_internal_slow_matches_[i]->compute(SemiAlignTypeFactory::getSuffixPtr(), suffixs);
    suffix_prsms_.push_back(suffixs);
    PrSMPtrVec internals;
    suffix_internal_slow_matches_[i]->compute(SemiAlignTypeFactory::getInternalPtr(), internals);
    internal_prsms_.push_back(internals);
  }
}

PrSMPtrVec PtmSlowFilter::getPrSMs(int nshift, SemiAlignTypePtr type){
    PrSMPtrVec matches;
  if (type == SemiAlignTypeFactory::getCompletePtr()) {
    for (unsigned int i = 0; i < complete_prsms_.size(); i++) {
      if (complete_prsms_[i][nshift] != nullptr) {
        matches.push_back(complete_prsms_[i][nshift]);
      }
    }
  }
  else if (type == SemiAlignTypeFactory::getPrefixPtr()) {
    for (unsigned int i = 0; i < prefix_prsms_.size(); i++) {
      if (prefix_prsms_[i][nshift] != nullptr) {
        matches.push_back(prefix_prsms_[i][nshift]);
      }
    }
  }
  else if (type == SemiAlignTypeFactory::getSuffixPtr()) {
    for (unsigned int i = 0; i < suffix_prsms_.size(); i++) {
      if (suffix_prsms_[i][nshift] != nullptr) {
        matches.push_back(suffix_prsms_[i][nshift]);
      }
    }
  }
  else if (type == SemiAlignTypeFactory::getInternalPtr()) {
    for (unsigned int i = 0; i < internal_prsms_.size(); i++) {
      if (internal_prsms_[i][nshift] != nullptr) {
        matches.push_back(internal_prsms_[i][nshift]);
      }
    }
  }
    return matches;
}

} /* namespace prot */
