//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#include "base/logger.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "poisson/poisson_comp_pvalue.hpp"

namespace prot {

PoissonCompPValue::PoissonCompPValue(ProteoformPtrVec &raw_forms, 
                                     ProteoformPtrVec &prot_mod_forms,
                                     PoissonMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  ResFreqPtrVec residue_freqs 
      = compResidueFreq(mng_ptr_->prsm_para_ptr_->getFixModResiduePtrVec(), raw_forms); 
  test_num_ptr_ = CountTestNumPtr(new CountTestNum(raw_forms, prot_mod_forms,
                                                   residue_freqs,
                                                   mng_ptr->convert_ratio_,
                                                   mng_ptr->max_prec_mass_,
                                                   mng_ptr->max_ptm_mass_));

  residue_avg_len_ = computeAvgLength(residue_freqs, 1.0);
  LOG_DEBUG("test number initialized")
}

double PoissonCompPValue::compRandMatchProb(double prec_mass, bool is_strict) {
  double avg_mass = prec_mass / 2.0;
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  // error tolerance for N-term + error tolerance for C -term
  double coverage;
  if (is_strict) {
    coverage = tole_ptr->compStrictErrorTole(avg_mass) * 2 + tole_ptr->compStrictErrorTole(avg_mass) * 2;
  }
  else {
    coverage = tole_ptr->compStrictErrorTole(avg_mass) * 2 + tole_ptr->compRelaxErrorTole(prec_mass, avg_mass) * 2;
  }
  double prob = coverage / residue_avg_len_;
  LOG_DEBUG(" Random match prob " << prob);
  return prob;
}

double PoissonCompPValue::compConditionProb(double rand_match_prob, ExtendMsPtr extend_ms, PrsmPtr prsm) {
  double f = extend_ms->size();
  double x = rand_match_prob;
  // use peak number is the same to the formula in prosightpc
  double n = prsm->getMatchPeakNum();
  LOG_DEBUG(" f " << f << " x " << x << " n " << n);
  n = n - prsm->getProteoformPtr()->getUnexpectedChangeNum() * mng_ptr_->shift_penalty;
  if (f <= 0.0 || x <= 0.0 || n <= 0.0) {
    return 1;
  }
  double mu = x*f;
  if (n <= mu) {
    return 1;
  }
  // approximation of poisson tail distribution
  double log_value = (n-mu) + n * std::log(mu/n) + std::log(n+1) - 0.5 * std::log(2*3.14 * n) - std::log(n+1-mu);
  double prob = std::exp(log_value);

  /*
  double xf = x*f;
  double sum = 0;
  double log_value = -xf;
  for (int i = 0; i < n; i++) {
    sum = sum + std::exp(log_value);
    log_value = log_value + std::log(xf) - std::log(i+1);
  }
  double prob = 1 - sum;
  */
  LOG_DEBUG(" f " << f << " x " << x  << " n " << n << " prob " << prob);
  return prob;
}

ExtremeValuePtrVec PoissonCompPValue::compExtremeValues(ExtendMsPtr extend_ms, 
                                                        PrsmPtrVec &prsms,
                                                        bool is_strict) {
  double prec_mass = extend_ms->getHeaderPtr()->getPrecMonoMassMinusWater();
  double tolerance = extend_ms->getHeaderPtr()->getErrorTolerance();
  double rand_match_prob = compRandMatchProb(prec_mass, is_strict); 

  ExtremeValuePtrVec ev_probs; 
  for (unsigned int i = 0; i < prsms.size(); i++) {
    double cond_prob = compConditionProb(rand_match_prob, extend_ms, prsms[i]);
    //LOG_DEBUG("prsm " << i << " prsm size " << prsms.size());
    int unexpect_shift_num = prsms[i]->getProteoformPtr()->getUnexpectedChangeNum();
    SemiAlignTypePtr t = prsms[i]->getProteoformPtr()->getSemiAlignType();
    double cand_num = test_num_ptr_->compCandNum(t, unexpect_shift_num, 
                                                 prec_mass, tolerance);
    //LOG_DEBUG("candidate number " << cand_num);
    if (cand_num == 0.0) {
      LOG_WARN("Zero candidate number");
      cand_num = std::numeric_limits<double>::infinity();
    }
    if (t == SemiAlignTypeFactory::getCompletePtr() 
        || t == SemiAlignTypeFactory::getPrefixPtr()) {
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(cond_prob, cand_num, 1)));
    } else {
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(cond_prob, cand_num,1)));
    }
    //LOG_DEBUG("assignment complete");
  }
  return ev_probs;
}

void PoissonCompPValue::setPValueArray(ExtendMsPtr extend_ms_ptr, PrsmPtrVec &prsms) {
  ExtremeValuePtrVec extreme_values = compExtremeValues(extend_ms_ptr, prsms, false);
  for (unsigned int i = 0; i < prsms.size(); i++) {
    prsms[i]->setProbPtr(extreme_values[i]);
  }
  //LOG_DEBUG("Set value complete");
}

}
