#ifndef ZERO_PTM_FILTER_COMP_SHIFT_HPP_
#define ZERO_PTM_FILTER_COMP_SHIFT_HPP_

#include <cmath>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"

namespace prot {

#define PRECURSOR_MATCH_SCORE 10000

class CompShift {
 public:
  CompShift(const ProteoformPtrVec &proteo_ptrs, int scale,
            double max_proteoform_mass, ProtModPtrVec prot_mod_ptr_vec);

  ~CompShift();

  void compZeroPtmConvolution(const std::vector<std::pair<int,int>> &mass_errors, 
                              std::pair<int,int> &prec_mass_error, 
                              int comp_num, int pref_suff_num, int inte_num);

  void compOnePtmConvolution(const std::vector<std::pair<int,int>> &mass_errors, 
                             int comp_num, int pref_suff_num, int inte_num);

  std::vector<std::pair<int,int>> getTopCompProteoScores() {return top_comp_proteo_scores_;}
  std::vector<std::pair<int,int>> getTopPrefProteoScores() {return top_pref_proteo_scores_;}
  std::vector<std::pair<int,int>> getTopSuffProteoScores() {return top_suff_proteo_scores_;}
  std::vector<std::pair<int,int>> getTopInternalProteoScores() {return top_internal_proteo_scores_;}

 private:
  // scale factor
  int scale_;
  //bool acetylation_;
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

  std::vector<int> col_index_begins_;
  std::vector<int> col_index_ends_;
  std::vector<int> col_indexes_;

  std::vector<int> rev_col_index_begins_;
  std::vector<int> rev_col_index_ends_;
  std::vector<int> rev_col_indexes_;

  std::vector<std::pair<int,int>> top_comp_proteo_scores_;
  std::vector<std::pair<int,int>> top_pref_proteo_scores_;
  std::vector<std::pair<int,int>> top_suff_proteo_scores_;
  std::vector<std::pair<int,int>> top_internal_proteo_scores_;

  void updateColumnMatchNums(ProteoformPtr proteo_ptr, ProtModPtr acet_mod, 
                             std::vector<int> &col_match_nums);
  void initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs);
  void initIndexes(const ProteoformPtrVec &proteo_ptrs);
  void updateRevColumnMatchNums(ProteoformPtr proteo_ptr, ProtModPtr acet_mod, 
                                std::vector<int> &col_match_nums);
  void initRevIndexes(const ProteoformPtrVec &proteo_ptrs);

  void compScores(const std::vector<std::pair<int,int>> &mass_errors,
                  std::vector<short> &scores, std::vector<short> &rev_scores);

  void findTopScores(std::vector<short> &scores, std::vector<short> &rev_scores, 
                     int comp_num, int pref_suff_num, int inte_num);
};

typedef std::shared_ptr<CompShift> CompShiftPtr;

} /* namespace prot */

#endif /* COMP_SHIFT_HPP_ */
