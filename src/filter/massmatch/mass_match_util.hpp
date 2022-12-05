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

#ifndef TOPPIC_FILTER_MASS_MATCH_MASS_MATCH_UTIL_HPP_
#define TOPPIC_FILTER_MASS_MATCH_MASS_MATCH_UTIL_HPP_

#include "seq/prot_candidate.hpp"
#include "filter/massmatch/mass_match.hpp"

namespace toppic {

namespace mass_match_util {

ProtCandidatePtrVec findZeroShiftTopProteins(std::vector<short> &scores,
                                             std::vector<short> &rev_scores,
                                             MassMatchPtr index_ptr,
                                             MassMatchPtr rev_index_ptr,
                                             double prec_minus_water_mass,
                                             double prec_error_tole, 
                                             int threshold, int num);

// the simple version is faster, but less accurate
ProtCandidatePtrVec simpleFindZeroShiftTopProteins(std::vector<short> &scores,
                                                   std::vector<short> &rev_scores,
                                                   MassMatchPtr index_ptr,
                                                   MassMatchPtr rev_index_ptr,
                                                   double prec_minus_water_mass,
                                                   double prec_error_tole, 
                                                   int threshold, int top_num);

ProtCandidatePtrVec findOneShiftTopProteins(std::vector<short> &scores,
                                            std::vector<short> &rev_scores,
                                            MassMatchPtr index_ptr,
                                            MassMatchPtr rev_index_ptr,
                                            double prec_minus_water_mass,
                                            double prec_error_tole,
                                            double min_shift, double max_shift,  
                                            int term_shift_num, 
                                            int threshold, int top_num); 

ProtCandidatePtrVec findVarPtmTopProteins(std::vector<short> &scores,
                                          std::vector<short> &rev_scores,
                                          MassMatchPtr index_ptr,
                                          MassMatchPtr rev_index_ptr,
                                          double prec_minus_water_mass,
                                          double prec_error_tole, 
                                          std::vector<double> &ptm_shifts,
                                          int threshold, int top_num); 

ProtCandidatePtrVec findDiagTopProteins(std::vector<short> &scores, 
                                        std::vector<int> &proteo_row_begins,
                                        std::vector<int> &proteo_row_ends,
                                        int threshold, int top_num);

}  // namespace mass_match_util

}  // namespace toppic

#endif 
