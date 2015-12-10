#ifndef ZERO_PTM_FILTER_ZERO_PTM_COMP_SHIFT_HPP_
#define ZERO_PTM_FILTER_ZERO_PTM_COMP_SHIFT_HPP_

#include <cmath>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "zeroptmfilter/zero_ptm_filter_mng.hpp"

namespace prot {

#define PRECURSOR_MATCH_SCORE 10000

class ZeroPtmCompShift {
 public:
  ZeroPtmCompShift(const ProteoformPtrVec &proteo_ptrs, ZeroPtmFilterMngPtr mng_ptr);

  ~ZeroPtmCompShift();

  void compConvolution(const std::vector<std::pair<int,int>> &mass_errors, 
                       std::pair<int,int> &prec_mass_error, ZeroPtmFilterMngPtr mng_ptr);

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
  int* proteo_row_begins_;
  int* proteo_row_ends_;
  ProtModPtr* acet_mods_;
  // the proteoform id of each row
  int* row_proteo_ids_;

  int* col_index_begins_;
  int* col_index_ends_;
  int* col_indexes_;

  int* rev_col_index_begins_;
  int* rev_col_index_ends_;
  int* rev_col_indexes_;

  std::vector<std::pair<int,int>> top_comp_proteo_scores_;
  std::vector<std::pair<int,int>> top_pref_proteo_scores_;
  std::vector<std::pair<int,int>> top_suff_proteo_scores_;
  std::vector<std::pair<int,int>> top_internal_proteo_scores_;

  void updateColumnMatchNums(ProteoformPtr proteo_ptr, ProtModPtr acet_mod, int* col_match_nums);
  void initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs);
  void initIndexes(const ProteoformPtrVec &proteo_ptrs);
  void updateRevColumnMatchNums(ProteoformPtr proteo_ptr, ProtModPtr acet_mod, int* col_match_nums);
  void initRevIndexes(const ProteoformPtrVec &proteo_ptrs);
  void compShiftScores(short* scores, short* rev_scores, ZeroPtmFilterMngPtr mng_ptr);
};

typedef std::shared_ptr<ZeroPtmCompShift> ZeroPtmCompShiftPtr;

} /* namespace prot */

#endif /* ZERO_PTM_COMP_SHIFT_HPP_ */
