#include "base/logger.hpp"
#include "tdgf/comp_pvalue_lookup_table.hpp"

namespace prot {

CompPValueLookupTable::CompPValueLookupTable(const ProteoformPtrVec &raw_proteo_ptrs, 
                                             const ProteoformPtrVec &mod_proteo_ptrs,
                                             const ResFreqPtrVec &residue_ptrs,
                                             TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  initTable();
  LOG_DEBUG("table initialized")
                        
  test_num_ptr_ = CountTestNumPtr(new CountTestNum(raw_proteo_ptrs, mod_proteo_ptrs,
                                                   residue_ptrs, mng_ptr->convert_ratio_,
                                                   mng_ptr->max_prec_mass_,
                                                   mng_ptr->max_ptm_mass_));
  LOG_DEBUG("test number initialized")
}

void CompPValueLookupTable::initTable() {
  // add init table
}

double CompPValueLookupTable::compProb(int ppo, double prec_mass, int peak_num, int match_frag_num, 
                                      int unexpected_shift_num) {
  // add implementation.
  return 1.0;
}

/* set alignment */
void CompPValueLookupTable::process(DeconvMsPtr deconv_ms_ptr, PrsmPtrVec &prsm_ptrs) {
  int ppo = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
  int peak_num = deconv_ms_ptr->size();
  double tolerance = deconv_ms_ptr->getHeaderPtr()->getErrorTolerance();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    double refine_prec_mass = prsm_ptrs[i]->getAdjustedPrecMass();
    int match_frag_num = prsm_ptrs[i]->getMatchFragNum();
    int unexpected_shift_num = prsm_ptrs[i]->getProteoformPtr()->getUnexpectedChangeNum();

    double prot_prob = compProb(ppo, refine_prec_mass, peak_num, match_frag_num, unexpected_shift_num);

    SemiAlignTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()->getSemiAlignType();

    double cand_num = test_num_ptr_->compCandNum(type_ptr, unexpected_shift_num, 
                                                 refine_prec_mass, tolerance);

    ExtremeValuePtr prob_ptr(new ExtremeValue(prot_prob, cand_num, 1));

    prsm_ptrs[i]->setProbPtr(prob_ptr);

  }
}

}
