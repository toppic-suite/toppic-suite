#ifndef PROT_COMP_PVALUE_LOOKUP_TABLE_HPP_
#define PROT_COMP_PVALUE_LOOKUP_TABLE_HPP_

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "tdgf/count_test_num.hpp"

namespace prot {

class CompPValueLookupTable {
 public:
  CompPValueLookupTable(TdgfMngPtr mng_ptr);

  void process(DeconvMsPtr deconv_ms_ptr, PrsmPtrVec &prsm_ptrs);

 private:

  void initTable();

  TdgfMngPtr mng_ptr_;
  CountTestNumPtr test_num_ptr_;
  std::ifstream input_;

  double ptm0_[41][20];
  double ptm1_[41][20];
  double ptm2_[41][20];

  double compProb(int match_frag_num, int unexpected_shift_num);
};

typedef std::shared_ptr<CompPValueLookupTable> CompPValueLookupTablePtr;

int getPeakIndex(int i);

int getFragIndex(int i);

}

#endif
