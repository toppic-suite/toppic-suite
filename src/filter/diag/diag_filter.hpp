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

#ifndef TOPPIC_FILTER_DIAG_DIAG_FILTER_HPP_
#define TOPPIC_FILTER_DIAG_DIAG_FILTER_HPP_

#include <memory>
#include <utility>
#include <vector>

#include "filter/massmatch/mass_match.hpp"
#include "filter/mng/diag_filter_mng.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

// A spectrum set converted once for the filter: the sorted integer masses
// with error windows that the index scores (prm_ms_util::getIntMassErrorList
// of the PRM spectra) and what the candidates need from the set (the header
// of its first spectrum and the number of spectra). It does not depend on
// the database block, so the block threads share one conversion.
class DiagFilterSpectrum {
 public:
  DiagFilterSpectrum(MsHeaderPtr header_ptr, int spectrum_num,
                     std::vector<std::pair<int, int>> mass_errors)
      : header_ptr_(std::move(header_ptr)),
        spectrum_num_(spectrum_num),
        mass_errors_(std::move(mass_errors)) {}

  const MsHeaderPtr& getHeaderPtr() const { return header_ptr_; }

  int getSpectrumNum() const { return spectrum_num_; }

  const std::vector<std::pair<int, int>>& getMassErrors() const {
    return mass_errors_;
  }

 private:
  MsHeaderPtr header_ptr_;
  int spectrum_num_;
  std::vector<std::pair<int, int>> mass_errors_;
};

using DiagFilterSpectrumPtr = std::shared_ptr<DiagFilterSpectrum>;
using DiagFilterSpectrumPtrVec = std::vector<DiagFilterSpectrumPtr>;

class DiagFilter {
 public:
  DiagFilter(const ProteoformPtrVec& proteo_ptrs,
             const DiagFilterMngPtr& mng_ptr, std::string block_str);

  SimplePrsmPtrVec getBestMatch(const PrmMsPtrVec& ms_ptr_vec);

  // The same candidates for a spectrum set converted by convert().
  SimplePrsmPtrVec getBestMatch(const DiagFilterSpectrum& spectrum);

  // Converts a spectrum set to what getBestMatch needs from it.
  static DiagFilterSpectrumPtr convert(const PrmMsPtrVec& ms_ptr_vec,
                                       const DiagFilterMngPtr& mng_ptr);

 private:
  DiagFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  MassMatchPtr index_ptr_;

  SimplePrsmPtrVec compute(const DiagFilterSpectrum& spectrum);
};

using DiagFilterPtr = std::shared_ptr<DiagFilter>;
} /* namespace toppic */

#endif
