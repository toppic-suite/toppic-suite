// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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
#include "base/neutral_loss_base.hpp"
#include "base/proteoform.hpp"
#include "spec/theo_peak_factory.hpp"

namespace prot {

TheoPeakPtrVec TheoPeakFactory::geneTheoPeak(BpSpecPtr bp_spec_ptr, ActivationPtr activation_ptr, 
                                             NeutralLossPtr neutral_loss_ptr,
                                             double n_term_shift, double c_term_shift,
                                             int bgn, int end, double min_mass, double max_mass){
  TheoPeakPtrVec theo_peaks;
  BreakPointPtrVec bps = bp_spec_ptr->getBreakPointPtrVec();
  IonTypePtr n_ion_type_ptr = activation_ptr->getNIonTypePtr();
  int charge = 0;
  for(int i =bgn; i<=end;i++){
    double n_mass = bps[i]->getNTermMass(n_ion_type_ptr);
    double new_mass = n_mass + n_term_shift;
    if(new_mass >= min_mass && new_mass <= max_mass){
      IonPtr ion = IonPtr(new Ion(charge,i,i,n_ion_type_ptr, neutral_loss_ptr));
      TheoPeakPtr theo_peak 
          = TheoPeakPtr(new TheoPeak(ion,n_mass,n_term_shift));
      theo_peaks.push_back(theo_peak);
    }
  }

  IonTypePtr c_ion_type_ptr = activation_ptr->getCIonTypePtr();
  for(int i =bgn;i<=end;i++){
    double c_mass = bps[i]->getCTermMass(c_ion_type_ptr);
    double new_mass = c_mass + c_term_shift;
    if(new_mass >= min_mass && new_mass <= max_mass){
      IonPtr ion 
          = IonPtr(new Ion(charge,i,bps.size()-i-1,c_ion_type_ptr,neutral_loss_ptr));
      theo_peaks.push_back(TheoPeakPtr(new TheoPeak(ion,c_mass,c_term_shift)));
    }
  }
  std::sort(theo_peaks.begin(),theo_peaks.end(),TheoPeak::cmpPosInc);
  return theo_peaks;
}

TheoPeakPtrVec TheoPeakFactory::geneProteoformTheoPeak(ProteoformPtr proteoform_ptr, 
                                                       ActivationPtr activation_ptr,
                                                       double min_mass) {
  BpSpecPtr bp_ptr = proteoform_ptr->getBpSpecPtr();

  TheoPeakPtrVec all_peaks;
  SegmentPtrVec segments = proteoform_ptr->getSegmentPtrVec();
  for (size_t i = 0; i < segments.size(); i++) {
    NeutralLossPtr neutral_loss_ptr 
        = NeutralLossBase::getNeutralLossPtr_NONE();
    double max_mass = proteoform_ptr->getResSeqPtr()->getSeqMass() 
        + segments[i]->getPepNTermShift() + segments[i]->getPepCTermShift() - min_mass; 
    TheoPeakPtrVec  peaks = geneTheoPeak(bp_ptr, 
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
