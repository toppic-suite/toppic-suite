#ifndef PROT_POSSION_COMP_PVALUE_HPP_
#define PROT_POSSION_COMP_PVALUE_HPP_

#include "tdgf/count_test_num.hpp"
#include "possion/possion_mng.hpp"

namespace prot {

class PossionCompPValue {
 public:
  PossionCompPValue(ProteoformPtrVec &raw_forms, 
                    ProteoformPtrVec &prot_mod_forms,
                    PossionMngPtr mng_ptr);

  ExtremeValuePtrVec compExtremeValues(PrmMsPtr ms_six, 
                                       PrsmPtrVec &prsms, bool strict);
  void setPValueArray(PrmMsPtr prm_ms_ptr, PrsmPtrVec &prsms);

 private:
  PossionMngPtr mng_ptr_;
  CountTestNumPtr test_num_ptr_;
  double residue_avg_len_;

  double compRandMatchProb(double prec_mass, bool is_strict);
  double compConditionProb(double rand_match_prob, PrmMsPtr ms_six, PrsmPtr prsm);
};

typedef std::shared_ptr<PossionCompPValue> PossionCompPValuePtr;

}

#endif
