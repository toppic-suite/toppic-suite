/*
 * theo_peak.cpp
 *
 *  Created on: Nov 28, 2013
 *      Author: xunlikun
 */

#include <algorithm>
#include "base/proteoform.hpp"
#include "spec/theo_peak.hpp"

namespace prot {

TheoPeak::TheoPeak(const IonPtr &ion_ptr,double unmod_mass,
                   double shift):
    Peak(unmod_mass + shift, 1.0) {
      ion_ptr_ = ion_ptr;
      unmod_mass_ = unmod_mass;
      shift_ = shift;
    }

TheoPeakPtrVec getTheoPeak(const BpSpecPtr &bp_spec,
                           const ActivationPtr &type, 
                           const NeutralLossPtr &neutral_loss_ptr,
                           double n_term_shift,
                           double c_term_shift,
                           int bgn,
                           int end,
                           double min_mass,
                           double max_mass){
  TheoPeakPtrVec theo_peaks;
  BreakPointPtrVec bps = bp_spec->getBreakPointPtrVec();
  IonTypePtr n_ion_type_ptr = type->getNIonTypePtr();
  for(int i =bgn;i<=end;i++){
    double n_mass = bps[i]->getNTermMass(n_ion_type_ptr);
    double new_mass = n_mass + n_term_shift;
    if(new_mass >= min_mass && new_mass <= max_mass){
      IonPtr ion = IonPtr(new Ion(0,i,i,n_ion_type_ptr, neutral_loss_ptr));
      TheoPeakPtr theo_peak 
          = TheoPeakPtr(new TheoPeak(ion,n_mass,n_term_shift));
      theo_peaks.push_back(theo_peak);
    }
  }

  IonTypePtr c_ion_type_ptr = type->getCIonTypePtr();
  for(int i =bgn;i<=end;i++){
    double c_mass = bps[i]->getCTermMass(c_ion_type_ptr);
    double new_mass = c_mass + c_term_shift;
    if(new_mass >= min_mass && new_mass <= max_mass){
      IonPtr ion 
          = IonPtr(new Ion(0,i,bps.size()-i-1,c_ion_type_ptr,neutral_loss_ptr));
      theo_peaks.push_back(TheoPeakPtr(new TheoPeak(ion,c_mass,c_term_shift)));
    }
  }
  std::sort(theo_peaks.begin(),theo_peaks.end(),theoPeakUp);
  return theo_peaks;
}

TheoPeakPtrVec getProteoformTheoPeak(const ProteoformPtr &proteoform_ptr, 
                                     const ActivationPtr &activation_ptr,
                                     double min_mass) {
  BpSpecPtr bp_ptr = proteoform_ptr->getBpSpecPtr();

  TheoPeakPtrVec all_peaks;
  SegmentPtrVec segments = proteoform_ptr->getSegmentPtrVec();
  for (unsigned int i = 0; i < segments.size(); i++) {
    NeutralLossPtr neutral_loss_ptr 
        = NeutralLossFactory::getNeutralLossPtr_NONE();
    double max_mass = proteoform_ptr->getResSeqPtr()->getSeqMass() 
        + segments[i]->getPepNTermShift() + segments[i]->getPepCTermShift() - min_mass; 
    TheoPeakPtrVec  peaks = getTheoPeak(bp_ptr, 
                                        activation_ptr, 
                                        neutral_loss_ptr,
                                        segments[i]->getPepNTermShift(),
                                        segments[i]->getPepCTermShift(), 
                                        segments[i]->getLeftBpPos(),
                                        segments[i]->getRightBpPos(),
                                        min_mass,
                                        max_mass);
    all_peaks.insert(all_peaks.end(), peaks.begin(), peaks.end());
  }
  return all_peaks;
}

} /* namespace prot */
