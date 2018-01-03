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


#ifndef PROT_BASE_BP_SPEC_HPP_
#define PROT_BASE_BP_SPEC_HPP_

#include <vector>

#include "base/residue_seq.hpp"
#include "base/break_point.hpp"

namespace prot {

// break point spectrum
class BpSpec {
 public:
  explicit BpSpec(const ResSeqPtr &res_seq_ptr);

  const BreakPointPtrVec& getBreakPointPtrVec() {return break_point_ptr_vec_;}

  BreakPointPtr getBreakPointPtr(int i) {return break_point_ptr_vec_[i];}

  // Get neutral ion masses for a specific ion type
  std::vector<double> getBreakPointMasses(IonTypePtr ion_type_ptr);

  std::vector<double> getPrmMasses();

  std::vector<double> getSrmMasses();

  // Get rounded scaled neutral ion masses
  std::vector<int> getScaledMass(double scale, IonTypePtr ion_type_ptr);

  std::vector<int> getScaledPrmMasses(double scale);

  std::vector<int> getScaledSrmMasses(double scale);

 private:
  BreakPointPtrVec break_point_ptr_vec_;

  void initBreakPoints(const ResSeqPtr &req_seq_ptr);
};

typedef std::shared_ptr<BpSpec> BpSpecPtr;
typedef std::vector<BpSpecPtr> BpSpecPtrVec;

}  // namespace prot

#endif
