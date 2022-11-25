//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <limits>

#include "common/base/mod_base.hpp"
#include "search/diag/diag_header_util.hpp"

namespace toppic {

namespace diag_header_util {

// get the header corresponding to the top left corner in the spectral grid
DiagHeaderPtr getTopLeftCornerHeader() {
  double shift = 0;
  // n_term strict; c_term nostrict; prot n_term match; prot c_term no_match
  // pep n_term no_match; pep c_term no_match
  return std::make_shared<DiagHeader>(shift, true, false, true, false, false, false);
}

DiagHeaderPtr getBottomRightCornerHeader(double seq_mass,
                                         double prec_mass) {
  double shift = prec_mass - seq_mass;
  // n term nostrict, c_term strict, prot n_term no match ; prot c_term match
  // pep n_term no match, pep c_term no match
  return std::make_shared<DiagHeader>(shift, false, true, false, true, false, false);
}

void addCornerDiagonals(DiagHeaderPtrVec &n_extend_header_ptrs,
                        DiagHeaderPtrVec &c_extend_header_ptrs,
                        double seq_mass, double prec_mass) {
  // get top-left corner header in spectral grid (shift is 0)
  DiagHeaderPtr top_left_corner_header_ptr = diag_header_util::getTopLeftCornerHeader();
  n_extend_header_ptrs.push_back(top_left_corner_header_ptr);

  // get bottom-right corner header in the spectral grid
  DiagHeaderPtr bottom_right_corner_header_ptr
      = getBottomRightCornerHeader(seq_mass, prec_mass);
  c_extend_header_ptrs.push_back(bottom_right_corner_header_ptr);
}

int findSimilarShiftPos(const std::vector<double> &shifts, double s) {
  int best_pos = -1;
  double best_diff = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < shifts.size(); i++) {
    if (std::abs(shifts[i] - s) < best_diff) {
      best_pos = i;
      best_diff = std::abs(shifts[i] - s);
    }
  }
  return best_pos;
}

bool isExistHeader(const DiagHeaderPtrVec &header_ptrs, double shift) {
  for (size_t i = 0; i < header_ptrs.size(); i++) {
    if (std::abs(header_ptrs[i]->getProtNTermShift()- shift) <= 0.01) {
      return true;
    }
  }
  return false;
}


DiagHeaderPtr geneDiagHeaderPtr(int match_bgn, int match_end, 
                                DiagHeaderPtr header_ptr) {

  DiagHeaderPtr new_header_ptr = std::make_shared<DiagHeader>(*header_ptr);
  new_header_ptr->setMatchFirstBpPos(match_bgn);
  new_header_ptr->setMatchLastBpPos(match_end);
  return new_header_ptr;
}

MassShiftPtrVec getDiagonalMassChanges(const DiagHeaderPtrVec &header_ptrs,
                                       int first_res_pos, int last_res_pos,
                                       AlterTypePtr type_ptr) {
  MassShiftPtrVec shift_list;
  ModPtr none_ptr = ModBase::getNoneModPtr();
  if (!header_ptrs[0]->isPepNTermMatch() && !header_ptrs[0]->isProtNTermMatch()) {
    AlterPtr alter_ptr
        = std::make_shared<Alter>(0, header_ptrs[0]->getMatchFirstBpPos()-first_res_pos,
                                  type_ptr, header_ptrs[0]->getPepNTermShift(), none_ptr);
    MassShiftPtr shift = std::make_shared<MassShift>(alter_ptr);
    shift_list.push_back(shift);
  }

  for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
    AlterPtr alter_ptr
        = std::make_shared<Alter>(header_ptrs[i]->getMatchLastBpPos() - first_res_pos,
                                  header_ptrs[i + 1]->getMatchFirstBpPos() - first_res_pos,
                                  type_ptr,
                                  header_ptrs[i + 1]->getProtNTermShift() 
                                  - header_ptrs[i]->getProtNTermShift(),
                                  none_ptr);
    MassShiftPtr shift = std::make_shared<MassShift>(alter_ptr);
    shift_list.push_back(shift);
  }

  DiagHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size() - 1];
  if (!last_header_ptr->isPepCTermMatch() && !last_header_ptr->isProtCTermMatch()) {
    AlterPtr alter_ptr
        = std::make_shared<Alter>(last_header_ptr->getMatchLastBpPos() - first_res_pos,
                                  (last_res_pos + 1) - first_res_pos,
                                  type_ptr, last_header_ptr->getPepCTermShift(),
                                  none_ptr);
    MassShiftPtr shift = std::make_shared<MassShift>(alter_ptr);
    shift_list.push_back(shift);
  }
  return shift_list;
}

MassShiftPtrVec getDiagonalMassChanges(const DiagHeaderPtrVec &header_ptrs,
                                       int first_res_pos, int last_res_pos,
                                       const AlterTypePtrVec & types) {
  MassShiftPtrVec shift_list;
  ModPtr none_ptr = ModBase::getNoneModPtr();
  if (!header_ptrs[0]->isPepNTermMatch() && !header_ptrs[0]->isProtNTermMatch()) {
    AlterPtr alter_ptr
        = std::make_shared<Alter>(0, header_ptrs[0]->getMatchFirstBpPos() - first_res_pos,
                                  types[0],
                                  header_ptrs[0]->getPepNTermShift(), none_ptr);

    MassShiftPtr shift = std::make_shared<MassShift>(alter_ptr);
    shift_list.push_back(shift);
  }

  for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
    AlterPtr alter_ptr
        = std::make_shared<Alter>(header_ptrs[i]->getMatchLastBpPos() - first_res_pos,
                                  header_ptrs[i + 1]->getMatchFirstBpPos() - first_res_pos,
                                  types[i + 1],
                                  header_ptrs[i + 1]->getProtNTermShift() 
                                  - header_ptrs[i]->getProtNTermShift(),
                                  none_ptr);
    MassShiftPtr shift = std::make_shared<MassShift>(alter_ptr);
    shift_list.push_back(shift);
  }

  DiagHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size() - 1];

  if (!last_header_ptr->isPepCTermMatch() && !last_header_ptr->isProtCTermMatch()) {
    AlterPtr alter_ptr
        = std::make_shared<Alter>(last_header_ptr->getMatchLastBpPos() - first_res_pos,
                                  (last_res_pos + 1) - first_res_pos,
                                  AlterType::UNEXPECTED,
                                  last_header_ptr->getPepCTermShift(),
                                  none_ptr);
    MassShiftPtr shift = std::make_shared<MassShift>(alter_ptr);
    shift_list.push_back(shift);
  }
  return shift_list;
}

}
}  // namespace toppic

