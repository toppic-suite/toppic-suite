// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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



#include "base/ion_type_base.hpp"
#include "spec/prm_peak_factory.hpp"

namespace prot {

PrmPeakPtr PrmPeakFactory::getZeroPeakPtr(int spec_id, double prec_mono_mass, 
                                          PeakTolerancePtr tole_ptr, double score) {
  //zero_peak
  DeconvPeakPtr zero_peak_ptr = DeconvPeakPtr(new DeconvPeak(-1,0,0,0));
  PrmPeakPtr prm_peak_ptr(new PrmPeak(spec_id, zero_peak_ptr,
                                      BasePeakType::ORIGINAL,0, score));
  // set tolerance 
  prm_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(0));
  prm_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(0));
  prm_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  return prm_peak_ptr;
}

PrmPeakPtr PrmPeakFactory::getPrecPeakPtr(int spec_id, double prec_mono_mass, 
                                          PeakTolerancePtr tole_ptr, double score) {
  //prec_peak
  double prec_peak_shift = IonTypeBase::getIonTypePtr_PREC()->getShift();
  double prec_peak_mass = prec_mono_mass - prec_peak_shift;
  DeconvPeakPtr prec_peak_ptr(new DeconvPeak(-1,prec_peak_mass, 0,0));
  PrmPeakPtr prm_peak_ptr = PrmPeakPtr(
      new PrmPeak(spec_id,prec_peak_ptr,
                  BasePeakType::ORIGINAL,prec_peak_mass, score));
  // set tolerance
  prm_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  prm_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  prm_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(0));
  return prm_peak_ptr;
}

} /* namespace prot */
