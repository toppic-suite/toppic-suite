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

#ifndef TOPPIC_FILTER_MASSMATCH_FILTER_SPEC_HPP_
#define TOPPIC_FILTER_MASSMATCH_FILTER_SPEC_HPP_

#include <memory>
#include <utility>
#include <vector>

#include "ms/spec/ms_header.hpp"

namespace toppic {

// The part of a spectrum set that a mass match filter uses: the integer mass
// and error lists of its prefix and suffix masses, the precursor mass error
// and the header of the first spectrum of the group. None of it depends on the
// database block, so a filter stage computes it once per spectrum set and
// shares it between the blocks instead of reading the msalign file and
// converting the spectra again for each of them.
class FilterSpec {
 public:
  FilterSpec(MsHeaderPtr header_ptr, int group_spec_num,
             std::vector<std::pair<int, int>> pref_mass_errors,
             std::vector<std::pair<int, int>> suff_mass_errors,
             std::pair<int, int> prec_minus_water_mass_error)
      : header_ptr_(std::move(header_ptr)),
        group_spec_num_(group_spec_num),
        pref_mass_errors_(std::move(pref_mass_errors)),
        suff_mass_errors_(std::move(suff_mass_errors)),
        prec_minus_water_mass_error_(prec_minus_water_mass_error) {}

  const MsHeaderPtr& getHeaderPtr() const { return header_ptr_; }

  int getGroupSpecNum() const { return group_spec_num_; }

  const std::vector<std::pair<int, int>>& getPrefMassErrors() const {
    return pref_mass_errors_;
  }

  const std::vector<std::pair<int, int>>& getSuffMassErrors() const {
    return suff_mass_errors_;
  }

  const std::pair<int, int>& getPrecMinusWaterMassError() const {
    return prec_minus_water_mass_error_;
  }

 private:
  MsHeaderPtr header_ptr_;
  int group_spec_num_;
  std::vector<std::pair<int, int>> pref_mass_errors_;
  std::vector<std::pair<int, int>> suff_mass_errors_;
  std::pair<int, int> prec_minus_water_mass_error_;
};

typedef std::shared_ptr<FilterSpec> FilterSpecPtr;
typedef std::vector<FilterSpecPtr> FilterSpecPtrVec;

}  // namespace toppic

#endif
