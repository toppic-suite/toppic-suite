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


#include <sstream>
#include <string>
#include <algorithm>

#include "seq/residue_seq.hpp"
#include "util/string_util.hpp"

namespace toppic {

ResidueSeq::ResidueSeq(const ResiduePtrVec &residues): residues_(residues) {
  /* get residue mass sum */
  residue_mass_sum_ = 0;
  for (size_t i = 0; i < residues_.size(); i++) {
    residue_mass_sum_ += residues_[i]->getMass();
  }
}

ResSeqPtr ResidueSeq::getSubResidueSeq(int bgn, int end) {
  if (end - bgn < 0) {
    return getEmptyResidueSeq();
  } else {
    ResiduePtrVec sub_residues;
    // from bgn to end,the sum of residues shoule be end - bgn + 1
    std::copy(residues_.begin() + bgn, residues_.begin() + end + 1,
              std::back_inserter(sub_residues) );
    return std::make_shared<ResidueSeq>(sub_residues);
  }
}

std::string ResidueSeq::toString() {
  std::stringstream s;
  for (size_t i = 0; i < residues_.size(); i++) {
    s << residues_[i]->toString();
  }
  s<< std::endl;
  return s.str();
}

std::string ResidueSeq::toAcidString() {
  std::stringstream s;
  for (size_t i = 0; i < residues_.size(); i++) {
    s << residues_[i]->getAminoAcidPtr()->getOneLetter();
  }
  return s.str();
}

ResSeqPtr ResidueSeq::getEmptyResidueSeq() {
  ResiduePtrVec residues;
  return std::make_shared<ResidueSeq>(residues);
}

}  // namespace toppic
