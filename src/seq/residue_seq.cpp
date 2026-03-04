//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#include "common/base/mod_base.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "seq/residue_seq.hpp"

namespace toppic {

ResidueSeq::ResidueSeq(const ResiduePtrVec &residues): 
    residues_(residues) {
  // get residue mass sum 
  residue_mass_sum_ = 0;
  for (size_t i = 0; i < residues_.size(); i++) {
    residue_mass_sum_ += residues_[i]->getMass();
  }
  n_mod_ptr_ = ModBase::getNTermNoneModPtr();
}

ResidueSeq::ResidueSeq(const ResiduePtrVec &residues, ModPtr n_mod_ptr): 
    residues_(residues), n_mod_ptr_(n_mod_ptr) {
  // get residue mass sum 
  residue_mass_sum_ = n_mod_ptr_->getShift();
  for (size_t i = 0; i < residues_.size(); i++) {
    residue_mass_sum_ += residues_[i]->getMass();
  }
}


ResSeqPtr ResidueSeq::getSubResidueSeq(int bgn, int end) {
  if (end - bgn < 0) {
    LOG_WARN("Empty sub residue sequence!");
    return getEmptyResidueSeq();
  } 
  if (end >= getLen()) {
    LOG_ERROR("The end posisiton is too large: " << end);
    exit(EXIT_FAILURE);
  }

  ResiduePtrVec sub_residues;
  // from bgn to end,the sum of residues shoule be end - bgn + 1
  std::copy(residues_.begin() + bgn, residues_.begin() + end + 1,
            std::back_inserter(sub_residues) );
  ModPtr n_mod_ptr_ = ModBase::getNTermNoneModPtr();
  if (bgn == 0) {
    n_mod_ptr_ = getNModPtr();
  }
  return std::make_shared<ResidueSeq>(sub_residues, n_mod_ptr_);
}

std::string ResidueSeq::toString() {
  std::stringstream s;
  if (!ModBase::isNTermNoneModPtr(n_mod_ptr_)) {
    s << "[" << n_mod_ptr_->getModResiduePtr()->getPtmPtr()->getAbbrName() << "]-"; 
  }
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

void ResidueSeq::setNModPtr(ModPtr n_mod_ptr) {
  n_mod_ptr_ = n_mod_ptr;
  residue_mass_sum_ = n_mod_ptr_->getShift();
  for (size_t i = 0; i < residues_.size(); i++) {
    residue_mass_sum_ += residues_[i]->getMass();
  }
}


}  // namespace toppic
