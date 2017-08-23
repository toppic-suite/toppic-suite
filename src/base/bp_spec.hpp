// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_BASE_BP_SPEC_HPP_
#define PROT_BASE_BP_SPEC_HPP_

#include <vector>

#include "base/residue_seq.hpp"
#include "base/break_point.hpp"

namespace prot {

// break point spectrum
class BpSpec {
 public:
  BpSpec() {}

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
