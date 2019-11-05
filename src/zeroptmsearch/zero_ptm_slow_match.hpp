//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_ZERO_PTM_SEARCH_ZERO_PTM_SLOW_MATCH_HPP_
#define TOPPIC_ZERO_PTM_SEARCH_ZERO_PTM_SLOW_MATCH_HPP_

#include "spec/deconv_peak.hpp"
#include "spec/theo_peak.hpp"
#include "zeroptmsearch/zero_ptm_search_mng.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"
#include "prsm/prsm.hpp"

namespace toppic {

class ZeroPtmSlowMatch;
typedef std::shared_ptr<ZeroPtmSlowMatch> ZpSlowMatchPtr;
typedef std::vector<ZpSlowMatchPtr> ZpSlowMatchPtrVec;

class ZeroPtmSlowMatch {
 public:
  ZeroPtmSlowMatch(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                   ZpFastMatchPtr fast_match_ptr,
                   ZeroPtmSearchMngPtr mng_ptr);

  ZeroPtmSlowMatch(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                   ProteoformPtr proteoform_ptr,
                   ZeroPtmSearchMngPtr mng_ptr);

  double getScore() {return score_;}

  PrsmPtr geneResult();

  static ZpSlowMatchPtrVec filter(const DeconvMsPtrVec &deconv_ms_ptr_vec,
                                  const ZpFastMatchPtrVec &fast_match_ptrs,
                                  ZeroPtmSearchMngPtr mng_ptr); 

 private:
  ZeroPtmSearchMngPtr mng_ptr_;
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  ZpFastMatchPtr fast_match_ptr_;
  ProteoformPtr proteoform_ptr_;

  double refine_prec_mass_;
  ExtendMsPtrVec refine_ms_ptr_vec_;

  double score_ = 0;

  void compScore (const ExtendMsPtrVec &refine_ms_ptr_vec);

  static bool cmpScoreDec(const ZpSlowMatchPtr &a, 
                          const ZpSlowMatchPtr &b) {
    return a->getScore() > b->getScore();
  }

  //double recal_ = 0;
  //bool isValid (double recal, double ppo);
};

}

#endif
