//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.



#include "base/ion_type_base.hpp"
#include "spec/prm_peak_factory.hpp"

namespace prot {

PrmPeakPtr PrmPeakFactory::getZeroPeakPtr(int spec_id, double prec_mono_mass,
                                          PeakTolerancePtr tole_ptr, double score) {
  // zero_peak
  DeconvPeakPtr zero_peak_ptr = std::make_shared<DeconvPeak>(-1, 0, 0, 0);
  PrmPeakPtr prm_peak_ptr = std::make_shared<PrmPeak>(spec_id, zero_peak_ptr,
                                                      BasePeakType::ORIGINAL, 0, score);
  // set tolerance
  prm_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(0));
  prm_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(0));
  prm_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  return prm_peak_ptr;
}

PrmPeakPtr PrmPeakFactory::getPrecPeakPtr(int spec_id, double prec_mono_mass,
                                          PeakTolerancePtr tole_ptr, double score) {
  // prec_peak
  double prec_peak_shift = IonTypeBase::getIonTypePtr_PREC()->getShift();
  double prec_peak_mass = prec_mono_mass - prec_peak_shift;
  DeconvPeakPtr prec_peak_ptr = std::make_shared<DeconvPeak>(-1, prec_peak_mass, 0, 0);
  PrmPeakPtr prm_peak_ptr
      = std::make_shared<PrmPeak>(spec_id, prec_peak_ptr,
                                  BasePeakType::ORIGINAL, prec_peak_mass, score);
  // set tolerance
  prm_peak_ptr->setStrictTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  prm_peak_ptr->setNStrictCRelacTolerance(tole_ptr->compStrictErrorTole(prec_mono_mass));
  prm_peak_ptr->setNRelaxCStrictTolerance(tole_ptr->compStrictErrorTole(0));
  return prm_peak_ptr;
}

} /* namespace prot */
