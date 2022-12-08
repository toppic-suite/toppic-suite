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

#ifndef TOPPIC_SEARCH_DIAG_DIAGONAL_UTIL_HPP_
#define TOPPIC_SEARCH_DIAG_DIAGONAL_UTIL_HPP_

#include "prsm/peak_ion_pair.hpp"
#include "search/diag/diag_header.hpp"

namespace toppic {

namespace diagonal_util {

double refinePrecursorAndHeaderShift(ProteoformPtr proteo_ptr,
                                     const ExtendMsPtrVec &ms_three_ptr_vec,
                                     DiagHeaderPtrVec &header_ptrs,
                                     double ppo, double min_mass,
                                     double refine_prec_step_width);
// for ps align
DiagHeaderPtrVec refineHeadersBgnEnd(ProteoformPtr proteo_ptr,
                                     const ExtendMsPtrVec &ms_three_ptr_vec,
                                     const DiagHeaderPtrVec& heade_ptrs,
                                     double min_mass);
// for ptm align
DiagHeaderPtrVec refineVarPtmHeadersBgnEnd(ProteoformPtr proteo_ptr,
                                           const ExtendMsPtrVec &ms_three_ptr_vec,
                                           const DiagHeaderPtrVec& header_ptrs,
                                           double min_mass); 
// for graph alignment
DiagHeaderPtrVec2D refineHeadersBgnEnd(ProteoformPtr proteo_ptr,
                                       const ExtendMsPtrVec &ms_three_ptr_vec,
                                       const DiagHeaderPtrVec2D& header_ptrs_2d,
                                       const DiagHeaderPtrVec& header_ptrs_1d,
                                       double min_mass);

int getNewBgn(const PeakIonPairPtrVec& pair_ptrs);

int getNewEnd(const PeakIonPairPtrVec& pair_ptrs);

}

} /* namespace toppic */

#endif /* DIAGONAL_HPP_ */
