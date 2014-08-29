#ifndef PROT_COMP_PVALUE_LOOKUP_TABLE_HPP_
#define PROT_COMP_PVALUE_LOOKUP_TABLE_HPP_

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "tdgf/count_test_num.hpp"

namespace prot {

class CompPValueLookupTable {
 public:
  CompPValueLookupTable(const ProteoformPtrVec &raw_proteo_ptrs, 
                        const ProteoformPtrVec &mod_proteo_ptrs,
                        const ResFreqPtrVec &residue_ptrs,
                        TdgfMngPtr mng_ptr);

  void process(DeconvMsPtr deconv_ms_ptr, PrsmPtrVec &prsm_ptrs);


 private:
  TdgfMngPtr mng_ptr_;
  CountTestNumPtr test_num_ptr_;
  
  void initTable();

  double compProb(int ppo, double prec_mass, int peak_num, 
                  int match_frag_num, int unexpected_shift_num);
};

typedef std::shared_ptr<CompPValueLookupTable> CompPValueLookupTablePtr;

}

#endif
