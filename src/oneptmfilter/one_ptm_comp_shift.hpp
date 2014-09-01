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

  std::vector<std::pair<int,int>> compConvolution(std::vector<int> &masses, 
                                                  std::vector<int> &errors, int num);

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

  void updateColumnMatchNums(ProteoformPtr proteo_ptr, int* col_match_nums);
  void initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs);
  void initIndexes(const ProteoformPtrVec &proteo_ptrs);
  void updateRevColumnMatchNums(ProteoformPtr proteo_ptr, int* col_match_nums);
  void initRevIndexes(const ProteoformPtrVec &proteo_ptrs);
  std::vector<std::pair<int,int>> getShiftScores(short* scores, short* rev_scores, int num);
};

typedef std::shared_ptr<OnePtmCompShift> OnePtmCompShiftPtr;

} /* namespace prot */

#endif /* ONE_PTM_COMP_SHIFT_HPP_ */
