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

#include "ms/spec/deconv_ms_util.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "ms/spec/deconv_peak.hpp"
#include "ms/spec/ms.hpp"
#include "ms/spec/ms_header.hpp"

namespace toppic {

namespace deconv_ms_util {

// Returns a deep copy of the deconvoluted spectra in which the precursor mass
// has been refined to new_prec_mass. Each spectrum gets a new header (cloned
// from the original but carrying new_prec_mass) and a fresh, independent copy
// of every deconvoluted peak, so the input vector is left unchanged.
DeconvMsPtrVec getRefineMsPtrVec(const DeconvMsPtrVec& deconv_ms_ptr_vec,
                                 double new_prec_mass) {
  DeconvMsPtrVec result_ptrs;
  for (size_t m = 0; m < deconv_ms_ptr_vec.size(); m++) {
    DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[m];
    MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
    MsHeaderPtr header_ptr =
        MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mass);
    std::vector<DeconvPeakPtr> peak_ptr_list;
    for (size_t p = 0; p < deconv_ms_ptr->size(); p++) {
      DeconvPeakPtr ori_peak_ptr = deconv_ms_ptr->getPeakPtr(p);
      // * is a dereference operator
      DeconvPeakPtr new_peak_ptr =
          std::make_shared<DeconvPeak>(*ori_peak_ptr.get());
      peak_ptr_list.push_back(new_peak_ptr);
    }
    DeconvMsPtr ms_ptr = std::make_shared<Ms<DeconvPeakPtr> >(
        header_ptr, std::move(peak_ptr_list));
    result_ptrs.push_back(ms_ptr);
  }
  return result_ptrs;
}

DeconvPeakPtrVec combineCloseMasses(const DeconvPeakPtrVec& peak_ptr_list,
                                    double ppo) {
  std::size_t n = peak_ptr_list.size();
  if (ppo <= 0 || n < 2) {
    return peak_ptr_list;
  }
  // the masses in increasing order, and each mass's position in that order
  std::vector<std::size_t> by_mass(n);
  for (std::size_t i = 0; i < n; i++) {
    by_mass[i] = i;
  }
  std::sort(by_mass.begin(), by_mass.end(),
            [&peak_ptr_list](std::size_t a, std::size_t b) {
              double ma = peak_ptr_list[a]->getMonoMass();
              double mb = peak_ptr_list[b]->getMonoMass();
              return (ma == mb) ? (a < b) : (ma < mb);
            });
  std::vector<std::size_t> rank(n);
  for (std::size_t k = 0; k < n; k++) {
    rank[by_mass[k]] = k;
  }
  // the masses in decreasing intensity, so that the most intense mass of a
  // group seeds it and the group is at most two tolerances wide
  std::vector<std::size_t> by_inte(by_mass.size());
  for (std::size_t i = 0; i < n; i++) {
    by_inte[i] = i;
  }
  std::sort(by_inte.begin(), by_inte.end(),
            [&peak_ptr_list](std::size_t a, std::size_t b) {
              double ia = peak_ptr_list[a]->getIntensity();
              double ib = peak_ptr_list[b]->getIntensity();
              return (ia == ib) ? (a < b) : (ia > ib);
            });

  const std::size_t none = n;
  std::vector<std::size_t> seed_of(n, none);
  for (std::size_t s : by_inte) {
    if (seed_of[s] != none) {
      continue;
    }
    seed_of[s] = s;
    double mass = peak_ptr_list[s]->getMonoMass();
    double tolerance = ppo * mass;
    for (std::size_t k = rank[s]; k > 0; k--) {
      std::size_t j = by_mass[k - 1];
      if (mass - peak_ptr_list[j]->getMonoMass() > tolerance) {
        break;
      }
      if (seed_of[j] == none) {
        seed_of[j] = s;
      }
    }
    for (std::size_t k = rank[s] + 1; k < n; k++) {
      std::size_t j = by_mass[k];
      if (peak_ptr_list[j]->getMonoMass() - mass > tolerance) {
        break;
      }
      if (seed_of[j] == none) {
        seed_of[j] = s;
      }
    }
  }

  std::vector<double> mass_sum(n, 0.0);
  std::vector<double> inte_sum(n, 0.0);
  std::vector<std::size_t> group_size(n, 0);
  for (std::size_t i = 0; i < n; i++) {
    std::size_t s = seed_of[i];
    double inte = peak_ptr_list[i]->getIntensity();
    mass_sum[s] += peak_ptr_list[i]->getMonoMass() * inte;
    inte_sum[s] += inte;
    group_size[s]++;
  }

  DeconvPeakPtrVec result;
  result.reserve(n);
  for (std::size_t i = 0; i < n; i++) {
    if (seed_of[i] != i) {
      continue;
    }
    if (group_size[i] == 1) {
      result.push_back(peak_ptr_list[i]);
      continue;
    }
    DeconvPeakPtr seed = peak_ptr_list[i];
    double mass = (inte_sum[i] > 0) ? (mass_sum[i] / inte_sum[i])
                                    : seed->getMonoMass();
    result.push_back(std::make_shared<DeconvPeak>(
        seed->getSpId(), seed->getPeakId(), mass, inte_sum[i],
        seed->getCharge(), seed->getScore()));
  }
  return result;
}

DeconvMsPtrVec combineCloseMasses(const DeconvMsPtrVec& deconv_ms_ptr_vec,
                                  double ppo) {
  if (ppo <= 0) {
    return deconv_ms_ptr_vec;
  }
  DeconvMsPtrVec result;
  result.reserve(deconv_ms_ptr_vec.size());
  for (const DeconvMsPtr& ms_ptr : deconv_ms_ptr_vec) {
    DeconvPeakPtrVec peaks = combineCloseMasses(ms_ptr->getPeakPtrVec(), ppo);
    result.push_back(std::make_shared<Ms<DeconvPeakPtr> >(
        ms_ptr->getMsHeaderPtr(), std::move(peaks)));
  }
  return result;
}

void keepTopPeaks(DeconvMsPtrVec& deconv_ms_ptr_vec, size_t peak_num) {
  for (size_t m = 0; m < deconv_ms_ptr_vec.size(); m++) {
    DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[m];
    std::vector<DeconvPeakPtr> peak_ptr_list;
    for (size_t p = 0; p < deconv_ms_ptr->size(); p++) {
      if (p >= peak_num) {
        break;
      }
      DeconvPeakPtr peak_ptr = deconv_ms_ptr->getPeakPtr(p);
      peak_ptr_list.push_back(peak_ptr);
    }
    deconv_ms_ptr->setPeakPtrVec(std::move(peak_ptr_list));
  }
}

}  // namespace deconv_ms_util

}  // namespace toppic
