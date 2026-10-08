// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "para/sp_para.hpp"

#include <memory>
#include <string>
#include <vector>

#include "common/base/activation_base.hpp"
#include "common/base/mass_constant.hpp"

namespace toppic {

SpPara::SpPara(const std::string& activation_name, double n_term_label_mass,
               double prec_ppm, double frag_ppm) {
  if (activation_name != "FILE") {
    activation_ptr_ = ActivationBase::getActivationPtrByName(activation_name);
  }
  n_term_label_mass_ = n_term_label_mass;
  prec_tolerance_ptr_ = std::make_shared<PeakTolerance>(prec_ppm * 0.000001);
  frag_tolerance_ptr_ = std::make_shared<PeakTolerance>(frag_ppm * 0.000001);

  // the set of offsets used to expand the monoisotopic mass list
  setFragMassExtension(true);

  setZeroShiftPrecOneDaError(true);
  // the one-shift stages absorb a 1 Da error into their unknown mass shift
  setOneShiftPrecOneDaError(false);
  setVarPtmPrecOneDaError(false);
  // multi-shift stages absorb a 1 Da error into an unknown mass shift
  multi_shift_search_prec_error_vec_ = {0};
}

namespace {

std::vector<double> precErrorVec(bool allow_one_da) {
  if (allow_one_da) {
    double IM = mass_constant::getIsotopeMass();
    return {0, -IM, IM};
  }
  return {0};
}

}  // namespace

void SpPara::setZeroShiftPrecOneDaError(bool allow) {
  zero_shift_search_prec_error_vec_ = precErrorVec(allow);
}

void SpPara::setOneShiftPrecOneDaError(bool allow) {
  one_shift_search_prec_error_vec_ = precErrorVec(allow);
}

void SpPara::setVarPtmPrecOneDaError(bool allow) {
  var_ptm_search_prec_error_vec_ = precErrorVec(allow);
}

void SpPara::setFragMassExtension(bool extend) {
  frag_ext_offsets_ = precErrorVec(extend);
}

const std::vector<std::string>& SpPara::getFilterMassListModeNames() {
  static const std::vector<std::string> names = {"merge", "fragment",
                                                 "single"};
  return names;
}

bool SpPara::isFilterMassListModeName(const std::string& name) {
  const std::vector<std::string>& names = getFilterMassListModeNames();
  for (size_t i = 0; i < names.size(); i++) {
    if (names[i] == name) {
      return true;
    }
  }
  return false;
}

bool SpPara::setFilterMassListMode(const std::string& name) {
  const std::vector<std::string>& names = getFilterMassListModeNames();
  for (size_t i = 0; i < names.size(); i++) {
    if (names[i] == name) {
      filter_mass_list_mode_ = static_cast<FilterMassListMode>(i);
      return true;
    }
  }
  return false;
}

}  // namespace toppic
