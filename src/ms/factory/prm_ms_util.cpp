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

#include "ms/factory/prm_ms_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "ms/factory/prm_peak_factory.hpp"
#include "ms/spec/base_peak_type.hpp"

namespace toppic {

namespace prm_ms_util {

namespace {

bool massErrorUp(const std::pair<int, int>& a, const std::pair<int, int>& b) {
  return a.first < b.first;
}

std::pair<int, int> getMassError(const PrmPeakPtr& peak_ptr, double scale,
                                 bool n_strict, bool c_strict) {
  int m = static_cast<int>(std::round(peak_ptr->getPosition() * scale));
  int e = 0;

  if (n_strict && c_strict) {
    e = static_cast<int>(std::ceil(peak_ptr->getStrictTolerance() * scale));
  } else if (n_strict && !c_strict) {
    e = static_cast<int>(
        std::ceil(peak_ptr->getNStrictCRelaxTolerance() * scale));
  } else if (!n_strict && c_strict) {
    e = static_cast<int>(
        std::ceil(peak_ptr->getNRelaxCStrictTolerance() * scale));
  }

  std::pair<int, int> mass_error(m, e);
  return mass_error;
}

}  // namespace

std::vector<std::pair<int, int> > getIntMassErrorList(
    const PrmMsPtrVec& prm_ms_ptr_vec, const PeakTolerancePtr& tole_ptr,
    double scale, bool n_strict, bool c_strict,
    SpPara::FilterMassListMode mode) {
  bool prec_n_strict = true;
  bool prec_c_strict = true;
  double default_score = 1.0;
  std::vector<std::pair<int, int> > mass_errors;
  for (size_t i = 0; i < prm_ms_ptr_vec.size(); i++) {
    const PrmPeakPtrVec& peak_ptrs = prm_ms_ptr_vec[i]->getPeakPtrVec();
    // kFragment (the default): one entry per integer mass and base type. Of
    // the peaks of one type with the same integer mass, the one with the
    // largest error is kept; an ORIGINAL and a REVERSED peak with the same
    // integer mass stay two entries. Such a pair is mostly an ion and the
    // PRM of its complementary ion (b and P - y), two matches of one prefix
    // mass of the proteoform when the precursor mass is accurate. As one
    // entry it scored one match, but two when a bin boundary fell between
    // the two peaks (toppic_test experiment 138).
    // kMerge: one entry per integer mass, whatever the type (all peaks in
    // slot 0). kSingle: one entry per peak; the peaks of one type with the
    // same integer mass are mostly one fragment mass at several charge
    // states, which TopFD lists separately (toppic_test experiment 142,
    // amendment 1).
    if (mode == SpPara::FilterMassListMode::kSingle) {
      for (size_t j = 0; j < peak_ptrs.size(); j++) {
        mass_errors.push_back(
            getMassError(peak_ptrs[j], scale, n_strict, c_strict));
      }
    } else {
      bool by_type = mode == SpPara::FilterMassListMode::kFragment;
      int last_mass[2] = {-1, -1};
      size_t last_index[2] = {0, 0};
      for (size_t j = 0; j < peak_ptrs.size(); j++) {
        std::pair<int, int> cur_m_e =
            getMassError(peak_ptrs[j], scale, n_strict, c_strict);
        int t = 0;
        if (by_type &&
            peak_ptrs[j]->getBaseTypePtr() == BasePeakType::REVERSED) {
          t = 1;
        }
        if (cur_m_e.first != last_mass[t]) {
          last_mass[t] = cur_m_e.first;
          last_index[t] = mass_errors.size();
          mass_errors.push_back(cur_m_e);
        } else if (cur_m_e.second > mass_errors[last_index[t]].second) {
          mass_errors[last_index[t]] = cur_m_e;
        }
      }
    }
    // add zero mass for each spectrum to increase the score for zero mass
    double prec_mass =
        prm_ms_ptr_vec[i]->getMsHeaderPtr()->getFirstPrecMonoMass();
    PrmPeakPtr zero_prm_ptr =
        prm_peak_factory::getZeroPeakPtr(i, prec_mass, tole_ptr, default_score);
    mass_errors.push_back(
        getMassError(zero_prm_ptr, scale, n_strict, c_strict));
    // add prec mass for each spectrum
    PrmPeakPtr prec_prm_ptr =
        prm_peak_factory::getPrecPeakPtr(i, prec_mass, tole_ptr, default_score);
    mass_errors.push_back(
        getMassError(prec_prm_ptr, scale, prec_n_strict, prec_c_strict));
  }

  std::sort(mass_errors.begin(), mass_errors.end(), massErrorUp);
  return mass_errors;
}

PrmPeakPtrVec getPrmPeakPtrs(const PrmMsPtrVec& prm_ms_ptr_vec,
                             const PeakTolerancePtr& tole_ptr) {
  PrmPeakPtrVec peak_list;
  for (size_t i = 0; i < prm_ms_ptr_vec.size(); i++) {
    for (size_t j = 0; j < prm_ms_ptr_vec[i]->size(); j++) {
      peak_list.push_back(prm_ms_ptr_vec[i]->getPeakPtr(j));
    }
  }
  // add zero
  double prec_mass =
      prm_ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
  // use spec_id = 0 and score = group_spec_num (size of prm_ms_ptr_vec)
  PrmPeakPtr zero_prm_ptr = prm_peak_factory::getZeroPeakPtr(
      0, prec_mass, tole_ptr, prm_ms_ptr_vec.size());
  peak_list.push_back(zero_prm_ptr);
  // add prec mass
  PrmPeakPtr prec_prm_ptr = prm_peak_factory::getPrecPeakPtr(
      0, prec_mass, tole_ptr, prm_ms_ptr_vec.size());
  peak_list.push_back(prec_prm_ptr);
  std::sort(peak_list.begin(), peak_list.end(), PrmPeak::cmpPosInc);
  for (size_t i = 0; i < peak_list.size(); i++) {
    peak_list[i]->setPeakId(i);
  }
  return peak_list;
}

}  // namespace prm_ms_util
}  // namespace toppic
