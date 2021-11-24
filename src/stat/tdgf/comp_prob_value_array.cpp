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

#include "stat/tdgf/comp_prob_value_array.hpp"

namespace toppic {

namespace comp_prob_value_array {

int getMaxScore(const PrsmPtrVec &prsm_ptrs) {
  int score = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getMatchFragNum() > score) {
      score = prsm_ptrs[i]->getMatchFragNum();
    }
  }
  return score;
}

int getMaxShift(PrsmPtrVec prsm_ptrs) {
  int shift = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getProteoformPtr()->getAlterNum(AlterType::UNEXPECTED) > shift) {
      shift = prsm_ptrs[i]->getProteoformPtr()->getAlterNum(AlterType::UNEXPECTED);
    }
  }
  return shift;
}

void compProbArray(CompProbValuePtr comp_prob_ptr, 
                   const ResFreqPtrVec &n_term_residue_ptrs, 
                   const PrmPeakPtrVec2D &peak_ptr_2d, 
                   const PrsmPtrVec &prsm_ptrs, 
                   bool strict, 
                   double prob_prec_mass,
                   PeakTolerancePtr tole_ptr, 
                   std::vector<double> &results) {
  int max_score = getMaxScore(prsm_ptrs);
  int max_shift = getMaxShift(prsm_ptrs);
  comp_prob_ptr->compute(n_term_residue_ptrs, peak_ptr_2d, max_score, max_shift, 
                         strict, prob_prec_mass, tole_ptr);
  results.clear();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int shift_num = prsm_ptrs[i]->getProteoformPtr()->getAlterNum(AlterType::UNEXPECTED);
    int score = prsm_ptrs[i]->getMatchFragNum();
    results.push_back(comp_prob_ptr->getCondProb(shift_num, score));
  }
}

}

}
