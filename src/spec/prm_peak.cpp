#include <algorithm>
#include <iostream>

#include "base/logger.hpp"
#include "base/support_peak_type_base.hpp"
#include "base/base_data.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

PrmPeak::PrmPeak(DeconvPeakPtr base_peak_ptr, int spec_id,
                 PrmBaseTypePtr base_type,
                 double mono_mass, double score):
    Peak(mono_mass, base_peak_ptr->getIntensity()),
    base_peak_ptr_(base_peak_ptr),
    spec_id_(spec_id),
    base_type_(base_type),
    mono_mass_(mono_mass),
    score_(score),
    strict_tolerance_(0),
    n_strict_c_relax_tolerance_(0),
    n_relax_c_strict_tolerance_(0) {
    }

void PrmPeak::addNghbEdge(DeconvPeakPtr deconv_peak_ptr,double offset,
                          SPTypePtr peak_type_ptr,double score){
  score_ +=score;
  SupportPeakPtr support_peak_ptr(
      new SupportPeak(deconv_peak_ptr,offset,score, peak_type_ptr));
  neighbor_list_.push_back(support_peak_ptr);
}

PrmBreakTypePtr PrmPeak::getBreakType() {
  PrmBreakTypePtr bt_ptr = PrmBreakType::NONE;
  for(size_t i=0; i<neighbor_list_.size(); i++){
    if(neighbor_list_[i]->getPeakTypePtr() == 
       SPTypeBase::getSPTypePtr_N_TERM()) {
      if(bt_ptr == PrmBreakType::NONE) {
        bt_ptr = PrmBreakType::N_TERM;
      }
      else if(bt_ptr == PrmBreakType::C_TERM){
        bt_ptr = PrmBreakType::BOTH;
      }
    }
    else{
      if(bt_ptr == PrmBreakType::NONE){
        bt_ptr = PrmBreakType::C_TERM;
      }
      else if(bt_ptr == PrmBreakType::N_TERM){
        bt_ptr = PrmBreakType::BOTH;
      }
    }
  }
  return bt_ptr;
}

} /* namespace prot */
