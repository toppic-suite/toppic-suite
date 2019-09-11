//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#ifndef TOPPIC_ONE_PTM_SEARCH_DIAGONAL_HEADER_UTIL_HPP_
#define TOPPIC_ONE_PTM_SEARCH_DIAGONAL_HEADER_UTIL_HPP_

#include <limits>
#include <cmath>
#include <vector>

#include "oneptmsearch/diagonal_header.hpp"

namespace toppic {
namespace DiagonalHeaderUtil {
// get the header corresponding to the top left corner in the spectral grid
DiagonalHeaderPtr getTopLeftCornerHeader();

DiagonalHeaderPtr getBottomRightCornerHeader(double seq_mass, double prec_mass);

void addCornerDiagonals(DiagonalHeaderPtrVec &n_extend_header_ptrs,
                        DiagonalHeaderPtrVec &c_extend_header_ptrs,
                        double seq_mass, double prec_mass);

int findSimilarShiftPos(const std::vector<double> &shifts, double s);

bool isExistHeader(const DiagonalHeaderPtrVec &header_ptrs, double shift);
}  // namespace DiagonalHeaderUtil
}  // namespace toppic

#endif
