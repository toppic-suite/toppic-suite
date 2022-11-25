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

#ifndef TOPPIC_SEARCH_DIAG_DIAG_HEADER_UTIL_HPP_
#define TOPPIC_SEARCH_DIAG_DIAG_HEADER_UTIL_HPP_

#include "search/diag/diag_header.hpp"

namespace toppic {

namespace diag_header_util {

// get the header corresponding to the top left corner in the spectral grid
DiagHeaderPtr getTopLeftCornerHeader();

DiagHeaderPtr getBottomRightCornerHeader(double seq_mass, double prec_mass);

void addCornerDiagonals(DiagHeaderPtrVec &n_extend_header_ptrs,
                        DiagHeaderPtrVec &c_extend_header_ptrs,
                        double seq_mass, double prec_mass);

int findSimilarShiftPos(const std::vector<double> &shifts, double s);

bool isExistHeader(const DiagHeaderPtrVec &header_ptrs, double shift);

// generate (clone) a new diagonal header with new match_bgn and match_end
DiagHeaderPtr geneDiagHeaderPtr(int match_bgn, int match_end,
                                DiagHeaderPtr diag_ptr);

MassShiftPtrVec getDiagonalMassChanges(const DiagHeaderPtrVec &diag_ptrs,
                                       int first_res_pos, int last_res_pos,
                                       AlterTypePtr type_ptr);

MassShiftPtrVec getDiagonalMassChanges(const DiagHeaderPtrVec &header_ptrs,
                                       int first_res_pos, int last_res_pos,
                                       const AlterTypePtrVec & type_ptrs);

}  // namespace diag_header_util

}  // namespace toppic

#endif
