#ifndef ONE_PTM_COMP_SHIFT_HPP_
#define ONE_PTM_COMP_SHIFT_HPP_

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "oneptmfilter/one_ptm_filter_mng.hpp"

namespace prot {

class OnePtmCompShift {
 public:
  OnePtmCompShift(const ProteoformPtrVec &proteo_ptrs, OnePtmFilterMngPtr mng_ptr);

  ~OnePtmCompShift();

  void compConvolution(const std::vector<std::pair<int,int>> &mass_errors, int num);

  std::vector<std::pair<int,int>> getTopCompProteoScores() {return top_comp_proteo_scores_;}
  std::vector<std::pair<int,int>> getTopPrefProteoScores() {return top_pref_proteo_scores_;}
  std::vector<std::pair<int,int>> getTopSuffProteoScores() {return top_suff_proteo_scores_;}
  std::vector<std::pair<int,int>> getTopInternalProteoScores() {return top_internal_proteo_scores_;}

 private:
  // scale factor
  int scale_;
  bool acetylation_;
  int proteo_num_;
  int col_num_;

  int row_num_;
  // the first row of each proteoform  
  int* proteo_row_begins_;
  int* proteo_row_ends_;
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

  void updateColumnMatchNums(ProteoformPtr proteo_ptr, int* col_match_nums);
  void initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs);
  void initIndexes(const ProteoformPtrVec &proteo_ptrs);
  void updateRevColumnMatchNums(ProteoformPtr proteo_ptr, int* col_match_nums);
  void initRevIndexes(const ProteoformPtrVec &proteo_ptrs);
  void compShiftScores(short* scores, short* rev_scores, int num);
};

typedef std::shared_ptr<OnePtmCompShift> OnePtmCompShiftPtr;

} /* namespace prot */

#endif /* ONE_PTM_COMP_SHIFT_HPP_ */
