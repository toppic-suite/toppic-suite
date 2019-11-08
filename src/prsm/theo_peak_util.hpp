//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#ifndef TOPPIC_PRSM_THEO_PEAK_UTIL_HPP_
#define TOPPIC_PRSM_THEO_PEAK_UTIL_HPP_

#include "common/base/activation.hpp"
#include "seq/bp_spec.hpp"
#include "seq/proteoform.hpp"
#include "ms/spec/theo_peak.hpp"

namespace toppic {

namespace theo_peak_util {

std::vector<double> getTheoMassVec(const TheoPeakPtrVec &theo_peak_list);

TheoPeakPtrVec geneTheoPeak(BpSpecPtr bp_spec_ptr, ActivationPtr activation_ptr,
                            NeutralLossPtr neutral_loss_ptr,
                            double n_term_shift, double c_term_shift,
                            int bgn, int end, double min_mass, double max_mass);

TheoPeakPtrVec geneProteoformTheoPeak(ProteoformPtr proteoform_ptr,
                                      ActivationPtr activation_ptr,
                                      double min_mass);

}  // namespace theo_peak_util

}  // namespace toppic

#endif

