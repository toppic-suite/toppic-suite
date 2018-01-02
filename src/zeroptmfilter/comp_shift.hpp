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


#ifndef ZERO_PTM_FILTER_COMP_SHIFT_HPP_
#define ZERO_PTM_FILTER_COMP_SHIFT_HPP_

#include <cmath>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "zeroptmfilter/filter_protein.hpp"

namespace prot {

#define PRECURSOR_MATCH_SCORE 10000

class CompShift {
 public:
  CompShift(const ProteoformPtrVec &proteo_ptrs, int scale,
            double max_proteoform_mass, ProtModPtrVec prot_mod_ptr_vec, 
            bool use_rev);

  ~CompShift() {}

  void compZeroPtmConvolution(const std::vector<std::pair<int,int>> &pref_mass_errors, 
                              const std::vector<std::pair<int,int>> &suff_mass_errors,
                              std::pair<int,int> &prec_minus_water_mass_error, 
                              int comp_num, int pref_suff_num, int inte_num);

  void compOnePtmConvolution(const std::vector<std::pair<int,int>> &pref_mass_errors, 
                             const std::vector<std::pair<int,int>> &suff_mass_errors,
                             int comp_num, int pref_suff_num, int inte_num, int shift_num);

  void compDiagConvolution(const std::vector<std::pair<int,int>> &mass_errors, 
                           int start, int top_num);

  FilterProteinPtrVec getTopCompProts() {return top_comp_prots_;}
  FilterProteinPtrVec getTopPrefProts() {return top_pref_prots_;}
  FilterProteinPtrVec getTopSuffProts() {return top_suff_prots_;}
  FilterProteinPtrVec getTopInternalProts() {return top_internal_prots_;}

  FilterProteinPtrVec getTopDiagProts() {return top_diag_prots_;}

 private:
  int scale_;
  ProtModPtrVec prot_mod_ptr_vec_;

  int proteo_num_;
  int col_num_;
  int row_num_;

  // the first row of each proteoform  
  std::vector<int> proteo_row_begins_;
  std::vector<int> proteo_row_ends_;
  ProtModPtrVec acet_mods_;
  // the proteoform id of each row
  std::vector<int> row_proteo_ids_;
  std::vector<double> n_trunc_shifts_;
  std::vector<double> c_trunc_shifts_;

  std::vector<int> col_index_begins_;
  std::vector<int> col_index_ends_;
  std::vector<int> col_indexes_;

  std::vector<int> rev_col_index_begins_;
  std::vector<int> rev_col_index_ends_;
  std::vector<int> rev_col_indexes_;

  FilterProteinPtrVec top_comp_prots_;
  FilterProteinPtrVec top_pref_prots_;
  FilterProteinPtrVec top_suff_prots_;
  FilterProteinPtrVec top_internal_prots_;

  FilterProteinPtrVec top_diag_prots_;

  void updateColumnMatchNums(ProteoformPtr proteo_ptr, ProtModPtr acet_mod, 
                             std::vector<int> &col_match_nums);
  void initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs);
  void initIndexes(const ProteoformPtrVec &proteo_ptrs);
  void updateRevColumnMatchNums(ProteoformPtr proteo_ptr, ProtModPtr acet_mod, 
                                std::vector<int> &col_match_nums);
  void initRevIndexes(const ProteoformPtrVec &proteo_ptrs);

  void compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                  std::vector<short> &scores);

  void compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                  int start, std::vector<short> &scores);

  void compRevScores(const std::vector<std::pair<int,int>> &suff_mass_errors,
                     std::vector<short> &rev_scores);

  void findTopScores(std::vector<short> &scores, std::vector<short> &rev_scores, 
                     double threshold, int comp_num, int pref_suff_num, int inte_num,
                     bool add_shifts, int shift_num);

  void findTopDiagScores(std::vector<short> &scores, int num);

  void addNTruncShifts(FilterProteinPtrVec &prot_ptrs, 
                       std::vector<short> &scores, int shift_num);

  void addCTruncShifts(FilterProteinPtrVec &prot_ptrs, 
                       std::vector<short> &scores, int shift_num);
};

typedef std::shared_ptr<CompShift> CompShiftPtr;

} /* namespace prot */

#endif /* COMP_SHIFT_HPP_ */
