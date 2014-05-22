#ifndef PROT_POSSION_COMP_PVALUE_HPP_
#define PROT_POSSION_COMP_PVALUE_HPP_

#include "tdgf/count_test_num.hpp"
#include "poisson/poisson_mng.hpp"

namespace prot {

class PoissonCompPValue {
 public:
  PoissonCompPValue(ProteoformPtrVec &raw_forms, 
                    ProteoformPtrVec &prot_mod_forms,
                    PoissonMngPtr mng_ptr);

  ExtremeValuePtrVec compExtremeValues(ExtendMsPtr extend_ms, 
                                       PrsmPtrVec &prsms, bool strict);
  void setPValueArray(ExtendMsPtr extend_ms_ptr, PrsmPtrVec &prsms);

 private:
  PoissonMngPtr mng_ptr_;
  CountTestNumPtr test_num_ptr_;
  double residue_avg_len_;

  double compRandMatchProb(double prec_mass, bool is_strict);
  double compConditionProb(double rand_match_prob, ExtendMsPtr extend_ms, PrsmPtr prsm);
};

typedef std::shared_ptr<PoissonCompPValue> PoissonCompPValuePtr;

}

#endif
