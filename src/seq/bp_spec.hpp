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

#ifndef TOPPIC_SEQ_BP_SPEC_HPP_
#define TOPPIC_SEQ_BP_SPEC_HPP_

#include <memory>
#include <vector>

#include "seq/break_point.hpp"
#include "seq/residue_seq.hpp"

namespace toppic {

// break point spectrum
class BpSpec {
 public:
  explicit BpSpec(const ResSeqPtr& res_seq_ptr);

  const BreakPointPtrVec& getBreakPointPtrVec() const {
    return break_point_ptr_vec_;
  }

  BreakPointPtr getBreakPointPtr(int i) const {
    return break_point_ptr_vec_[i];
  }

  // Get neutral ion masses for a specific ion type
  std::vector<double> getBreakPointMasses(const IonTypePtr& ion_type_ptr) const;

  std::vector<double> getPrmMasses() const;

  std::vector<double> getSrmMasses() const;

  // Get rounded scaled neutral ion masses
  std::vector<int> getScaledMass(double scale,
                                 const IonTypePtr& ion_type_ptr) const;

  std::vector<int> getScaledPrmMasses(double scale) const;

  std::vector<int> getScaledSrmMasses(double scale) const;

 private:
  BreakPointPtrVec break_point_ptr_vec_;

  void initBreakPoints(const ResSeqPtr& req_seq_ptr);
};

using BpSpecPtr = std::shared_ptr<BpSpec>;
using BpSpecPtrVec = std::vector<BpSpecPtr>;

}  // namespace toppic

#endif
