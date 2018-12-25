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


#include <vector>
#include <algorithm>

#include "seq/bp_spec.hpp"
#include "common/base/neutral_loss_base.hpp"
#include "prsm/theo_peak_util.hpp"

namespace toppic {

namespace theo_peak_util {

std::vector<double> getTheoMassVec(const TheoPeakPtrVec &theo_peak_list) {
  std::vector<double> masses;
  for (size_t i = 0; i < theo_peak_list.size(); i++) {
    masses.push_back(theo_peak_list[i]->getModMass());
  }
  return masses;
}

TheoPeakPtrVec geneTheoPeak(BpSpecPtr bp_spec_ptr, ActivationPtr activation_ptr,
                            NeutralLossPtr neutral_loss_ptr,
                            double n_term_shift, double c_term_shift,
                            int bgn, int end, double min_mass, double max_mass) {
  TheoPeakPtrVec theo_peaks;
  BreakPointPtrVec bps = bp_spec_ptr->getBreakPointPtrVec();
  IonTypePtr n_ion_type_ptr = activation_ptr->getNIonTypePtr();
  int charge = 0;
  for (int i = bgn; i <= end; i++) {
    double n_mass = bps[i]->getNTermMass(n_ion_type_ptr);
    double new_mass = n_mass + n_term_shift;
    if (new_mass >= min_mass && new_mass <= max_mass) {
      IonPtr ion = std::make_shared<Ion>(charge, i, i,
                                         n_ion_type_ptr, neutral_loss_ptr);
      TheoPeakPtr theo_peak = std::make_shared<TheoPeak>(ion, n_mass, n_term_shift);
      theo_peaks.push_back(theo_peak);
    }
  }

  IonTypePtr c_ion_type_ptr = activation_ptr->getCIonTypePtr();
  for (int i = bgn; i <= end; i++) {
    double c_mass = bps[i]->getCTermMass(c_ion_type_ptr);
    double new_mass = c_mass + c_term_shift;
    if (new_mass >= min_mass && new_mass <= max_mass) {
      IonPtr ion = std::make_shared<Ion>(charge, i, bps.size()-i-1,
                                         c_ion_type_ptr, neutral_loss_ptr);
      theo_peaks.push_back(std::make_shared<TheoPeak>(ion, c_mass, c_term_shift));
    }
  }
  std::sort(theo_peaks.begin(), theo_peaks.end(), TheoPeak::cmpPosInc);
  return theo_peaks;
}

TheoPeakPtrVec geneProteoformTheoPeak(ProteoformPtr proteoform_ptr,
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

}  // namespace theo_peak_util

}  // namespace toppic
