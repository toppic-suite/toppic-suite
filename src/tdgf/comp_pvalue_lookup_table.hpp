#ifndef PROT_COMP_PVALUE_LOOKUP_TABLE_HPP_
#define PROT_COMP_PVALUE_LOOKUP_TABLE_HPP_

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "tdgf/count_test_num.hpp"

namespace prot {

class CompPValueLookupTable {
 public:
  CompPValueLookupTable(TdgfMngPtr mng_ptr);

  bool inTable(const DeconvMsPtrVec &deconv_ms_ptr_vec, const PrsmPtrVec &prsm_ptrs);

  void process(const DeconvMsPtrVec &deconv_ms_ptr_vec, PrsmPtrVec &prsm_ptrs, double ppo);

 private:

  void initTable();

  TdgfMngPtr mng_ptr_;
  CountTestNumPtr test_num_ptr_;
  std::ifstream input_;

  double ptm0_[48][20];
  double ptm1_[48][20];
  double ptm2_[48][20];

  double compProb(int peak_num, int match_frag_num, int unexpected_shift_num);
};

typedef std::shared_ptr<CompPValueLookupTable> CompPValueLookupTablePtr;

int getPeakIndex(int i);

int getFragIndex(int i);

// get x1, x2, y1, y2
std::vector<int> getFourIndex(int peak_num, int frag_num);

int getPeakNumFromIndex(int idx);

}

#endif
