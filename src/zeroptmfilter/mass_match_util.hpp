//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef ZERO_PTM_FILTER_MASS_MATCH_UTIL_HPP_
#define ZERO_PTM_FILTER_MASS_MATCH_UTIL_HPP_

#include <vector>
#include <cmath>

#include "zeroptmfilter/filter_protein.hpp"
#include "zeroptmfilter/mass_match.hpp"

namespace prot {

namespace mass_match_util {

FilterProteinPtrVec findTopProteins(std::vector<short> &scores, 
                                    std::vector<int> &proteo_row_begins,
                                    std::vector<int> &proteo_row_ends,
                                    int threshold, int num);

FilterProteinPtrVec findTopProteins(std::vector<short> &scores, 
                                    std::vector<short> &rev_scores, 
                                    MassMatchPtr index_ptr,
                                    MassMatchPtr rev_index_ptr,
                                    double threshold, int num,
                                    bool add_shifts, int shift_num);
}  // namespace mass_match_util

}  // namespace prot

#endif 
