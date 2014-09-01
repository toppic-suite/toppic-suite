#ifndef PROT_COMP_SHIFT_HPP_
#define PROT_COMP_SHIFT_HPP_

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "diagfilter/diag_filter_mng.hpp"

namespace prot {

class CompShift {
 public:
  CompShift(const ProteoformPtrVec &proteo_ptrs, DiagFilterMngPtr mng_ptr);

  ~CompShift();

  std::vector<std::pair<int,int>> compConvolution(const std::vector<int> &masses,
                                                  int bgn_pos,int num);

  std::vector<std::pair<int,int>> compConvolution(const std::vector<int> &masses, 
                                                  const std::vector<int> &errors,
                                                  int bgn_pos,int num);

 private:
  // scale factor
  int scale_;
  int col_num_;

  int row_num_;
  // the first row of each proteoform  
  int* proteo_row_begins_;
  // the proteoform id of each row
  int* row_proteo_ids_;

  int* col_index_begins_;
  int* col_index_ends_;
  int* col_indexes_;

  void updateColumnMatchNums(ProteoformPtr proteo_ptr, int* col_match_nums, 
                             bool acetylation);
  void initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs, bool acetylation);
  void initIndexes(const ProteoformPtrVec &proteo_ptrs, bool acetylation);
  std::vector<std::pair<int,int>> getShiftScores(short* scores, int num);
};

typedef std::shared_ptr<CompShift> CompShiftPtr;

} /* namespace prot */

#endif /* PROT_COMP_SHIFT_HPP_ */
