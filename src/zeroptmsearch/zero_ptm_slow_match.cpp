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
#include "base/proteoform.hpp"
#include "base/proteoform_factory.hpp"
#include "base/activation.hpp"
#include "base/base_algo.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/theo_peak.hpp"
#include "spec/theo_peak_factory.hpp"
#include "spec/theo_peak_util.hpp"
#include "zeroptmsearch/zero_ptm_slow_match.hpp"

namespace prot {

ZeroPtmSlowMatch::ZeroPtmSlowMatch(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                   ZpFastMatchPtr fast_match_ptr,
                                   ZeroPtmSearchMngPtr mng_ptr): 
    mng_ptr_(mng_ptr),
    deconv_ms_ptr_vec_(deconv_ms_ptr_vec),
    fast_match_ptr_(fast_match_ptr) {

      proteoform_ptr_ = ProteoformFactory::geneSubProteoform(fast_match_ptr->getProteoformPtr(), 
                                                             fast_match_ptr->getBegin(), 
                                                             fast_match_ptr->getEnd());

      SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
      refine_prec_mass_ = proteoform_ptr_->getResSeqPtr()->getSeqMass();
      refine_ms_ptr_vec_ = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec_, 
                                                              sp_para_ptr, 
                                                              refine_prec_mass_);

      compScore(refine_ms_ptr_vec_);
    }

ZeroPtmSlowMatch::ZeroPtmSlowMatch(const DeconvMsPtrVec &deconv_ms_ptr_vec,                               
                                   ProteoformPtr proteoform_ptr,                                        
                                   ZeroPtmSearchMngPtr mng_ptr):
    mng_ptr_(mng_ptr),
    deconv_ms_ptr_vec_(deconv_ms_ptr_vec),
    proteoform_ptr_(proteoform_ptr) {
      SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
      refine_prec_mass_ = proteoform_ptr_->getResSeqPtr()->getSeqMass();
      refine_ms_ptr_vec_ = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec_, 
                                                              sp_para_ptr, 
                                                              refine_prec_mass_);
      compScore(refine_ms_ptr_vec_);
    }

// compute the average ppo
double compAvg(const std::vector<double> &ppos, double recal_ppo) {
  int cnt = 0;
  double sum = 0;
  for (size_t i = 0; i < ppos.size(); i++) {
    if (std::abs(ppos[i]) <= recal_ppo) {
      cnt++;
      sum += ppos[i];
    }
  }
  if (cnt == 0) {
    return 0;
  } else {
    return sum / cnt;
  }
}

// input is refineMsThree, result is score and recal and recalMass 
void ZeroPtmSlowMatch::compScore (const ExtendMsPtrVec &refine_ms_ptr_vec) {
  score_ = 0;
  double min_mass = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getMinMass();
  double ppo = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();
  for (size_t i = 0; i < refine_ms_ptr_vec.size(); i++) {
    ActivationPtr activation_ptr = refine_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr();
    TheoPeakPtrVec theo_peak_ptrs = TheoPeakFactory::geneProteoformTheoPeak(proteoform_ptr_, 
                                                                            activation_ptr, min_mass);

    std::vector<double> theo_masses = Theopeak_util::getTheoMassVec(theo_peak_ptrs);
    std::vector<double> ms_masses = ExtendMs::getExtendMassVec(refine_ms_ptr_vec[i]);
    score_ += base_algo::compNumMatchedTheoMasses(ms_masses, theo_masses, ppo);
  }
}

// get result 
PrsmPtr ZeroPtmSlowMatch::geneResult() {
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  return std::make_shared<Prsm>(proteoform_ptr_, deconv_ms_ptr_vec_,
                                refine_prec_mass_, sp_para_ptr);
}

ZpSlowMatchPtrVec ZeroPtmSlowMatch::filter(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                                           const ZpFastMatchPtrVec &fast_match_ptrs,
                                           ZeroPtmSearchMngPtr mng_ptr) {

  ZpSlowMatchPtrVec slow_matches;
  for (size_t i = 0; i < fast_match_ptrs.size(); i++) {
    ZpSlowMatchPtr slow_match
        = std::make_shared<ZeroPtmSlowMatch>(deconv_ms_ptr_vec, fast_match_ptrs[i], mng_ptr);
    slow_matches.push_back(slow_match);
  }
  // sort 
  std::sort(slow_matches.begin(), slow_matches.end(), ZeroPtmSlowMatch::cmpScoreDec);

  return slow_matches;
}

}
