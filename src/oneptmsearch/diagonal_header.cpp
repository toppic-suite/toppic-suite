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


#include <cmath>
#include <vector>

#include "base/prot_mod.hpp"
#include "base/mod_base.hpp"
#include "base/change.hpp"
#include "base/base_algo.hpp"
#include "oneptmsearch/diagonal_header.hpp"

namespace prot {

DiagonalHeaderPtr DiagonalHeader::clone() {
  DiagonalHeaderPtr cloned
      = std::make_shared<DiagonalHeader>(prot_N_term_shift_, n_strict_, c_strict_,
                                         prot_N_term_match_, prot_C_term_match_,
                                         pep_N_term_match_, pep_C_term_match_);
  cloned->setId(id_);
  cloned->setTruncFirstResPos(trunc_first_res_pos_);
  cloned->setMatchFirstBpPos(match_first_bp_pos_);
  cloned->setPepNTermShift(pep_N_term_shift_);
  cloned->setProtNTermShift(prot_N_term_shift_);
  cloned->setAlignPrefix(is_align_prefix_);
  cloned->setTruncLastResPos(trunc_last_res_pos_);
  cloned->setMatchLastBpPos(match_last_bp_pos_);
  cloned->setPepCTermShift(pep_C_term_shift_);
  cloned->setProtCTermShift(prot_C_term_shift_);
  cloned->setAlignSuffix(is_align_suffix_);
  return cloned;
}

DiagonalHeaderPtr geneDiagonalHeaderPtr(int bgn, int end,
                                        DiagonalHeaderPtr header_ptr) {
  DiagonalHeaderPtr new_header_ptr = header_ptr->clone();
  new_header_ptr->setMatchFirstBpPos(bgn);
  new_header_ptr->setMatchLastBpPos(end);
  return new_header_ptr;
}

void DiagonalHeader::initHeader(double c_shift, ProteoformPtr proteo_ptr,
                                double align_pref_suff_shift_thresh) {
  // set protein c term shift
  prot_C_term_shift_ = c_shift;
  std::vector<double> prm_masses = proteo_ptr->getBpSpecPtr()->getPrmMasses();

  trunc_first_res_pos_ = base_algo::getFirstResPos(prot_N_term_shift_,
                                                   prm_masses);

  trunc_last_res_pos_ = base_algo::getLastResPos(prot_C_term_shift_, prm_masses);
  pep_N_term_shift_ = prot_N_term_shift_ + prm_masses[trunc_first_res_pos_];
  pep_C_term_shift_ = prot_C_term_shift_ + prm_masses[prm_masses.size() - 1]
      - prm_masses[trunc_last_res_pos_ + 1];

  if (std::abs(prot_N_term_shift_) <= align_pref_suff_shift_thresh) {
    is_align_prefix_ = true;
  } else {
    is_align_prefix_ = false;
  }

  if (std::abs(prot_C_term_shift_) <= align_pref_suff_shift_thresh) {
    is_align_suffix_ = true;
  } else {
    is_align_suffix_ = false;
  }
}

MassShiftPtrVec getDiagonalMassChanges(const DiagonalHeaderPtrVec &header_ptrs,
                                       int first_res_pos, int last_res_pos,
                                       MassShiftTypePtr type_ptr) {
  MassShiftPtrVec shift_list;
  ModPtr none_ptr = ModBase::getNoneModPtr();
  if (!header_ptrs[0]->isPepNTermMatch() && !header_ptrs[0]->isProtNTermMatch()) {
    ChangePtr c
        = std::make_shared<Change>(0, header_ptrs[0]->getMatchFirstBpPos()-first_res_pos,
                                   type_ptr, header_ptrs[0]->getPepNTermShift(), none_ptr);
    MassShiftPtr shift
        = std::make_shared<MassShift>(c->getLeftBpPos(), c->getRightBpPos(), c->getTypePtr());
    shift->setChangePtr(c);
    shift_list.push_back(shift);
  }

  for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
    ChangePtr change_ptr
        = std::make_shared<Change>(header_ptrs[i]->getMatchLastBpPos() - first_res_pos,
                                   header_ptrs[i + 1]->getMatchFirstBpPos() - first_res_pos,
                                   type_ptr,
                                   header_ptrs[i + 1]->getProtNTermShift() - header_ptrs[i]->getProtNTermShift(),
                                   none_ptr);
    MassShiftPtr shift
        = std::make_shared<MassShift>(change_ptr->getLeftBpPos(),
                                      change_ptr->getRightBpPos(),
                                      change_ptr->getTypePtr());
    shift->setChangePtr(change_ptr);
    shift_list.push_back(shift);
  }

  DiagonalHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size() - 1];
  if (!last_header_ptr->isPepCTermMatch() && !last_header_ptr->isProtCTermMatch()) {
    ChangePtr change_ptr
        = std::make_shared<Change>(last_header_ptr->getMatchLastBpPos() - first_res_pos,
                                   (last_res_pos + 1) - first_res_pos,
                                   type_ptr, last_header_ptr->getPepCTermShift(),
                                   none_ptr);
    MassShiftPtr shift
        = std::make_shared<MassShift>(change_ptr->getLeftBpPos(),
                                      change_ptr->getRightBpPos(),
                                      change_ptr->getTypePtr());
    shift->setChangePtr(change_ptr);
    shift_list.push_back(shift);
  }
  return shift_list;
}

MassShiftPtrVec getDiagonalMassChanges(const DiagonalHeaderPtrVec &header_ptrs,
                                       int first_res_pos, int last_res_pos,
                                       const MassShiftTypePtrVec & types) {
  MassShiftPtrVec shift_list;
  ModPtr none_ptr = ModBase::getNoneModPtr();
  if (!header_ptrs[0]->isPepNTermMatch() && !header_ptrs[0]->isProtNTermMatch()) {
    ChangePtr change_ptr
        = std::make_shared<Change>(0, header_ptrs[0]->getMatchFirstBpPos() - first_res_pos,
                                   types[0],
                                   header_ptrs[0]->getPepNTermShift(), none_ptr);

    MassShiftPtr shift
        = std::make_shared<MassShift>(change_ptr->getLeftBpPos(),
                                      change_ptr->getRightBpPos(),
                                      change_ptr->getTypePtr());
    shift->setChangePtr(change_ptr);
    shift_list.push_back(shift);
  }

  for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
    ChangePtr change_ptr
        = std::make_shared<Change>(header_ptrs[i]->getMatchLastBpPos() - first_res_pos,
                                   header_ptrs[i + 1]->getMatchFirstBpPos() - first_res_pos,
                                   types[i + 1],
                                   header_ptrs[i + 1]->getProtNTermShift() - header_ptrs[i]->getProtNTermShift(),
                                   none_ptr);
    MassShiftPtr shift
        = std::make_shared<MassShift>(change_ptr->getLeftBpPos(),
                                      change_ptr->getRightBpPos(),
                                      change_ptr->getTypePtr());
    shift->setChangePtr(change_ptr);
    shift_list.push_back(shift);
  }

  DiagonalHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size() - 1];

  if (!last_header_ptr->isPepCTermMatch() && !last_header_ptr->isProtCTermMatch()) {
    ChangePtr change_ptr
        = std::make_shared<Change>(last_header_ptr->getMatchLastBpPos() - first_res_pos,
                                   (last_res_pos + 1) - first_res_pos,
                                   MassShiftType::UNEXPECTED,
                                   last_header_ptr->getPepCTermShift(),
                                   none_ptr);
    MassShiftPtr shift
        = std::make_shared<MassShift>(change_ptr->getLeftBpPos(),
                                      change_ptr->getRightBpPos(),
                                      change_ptr->getTypePtr());
    shift->setChangePtr(change_ptr);
    shift_list.push_back(shift);
  }
  return shift_list;
}

}  // namespace prot
