// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef ZERO_PTM_SLOW_MATCH_HPP_
#define ZERO_PTM_SLOW_MATCH_HPP_

#include "spec/deconv_peak.hpp"
#include "spec/theo_peak.hpp"
#include "zeroptmsearch/zero_ptm_search_mng.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class ZeroPtmSlowMatch;
typedef std::shared_ptr<ZeroPtmSlowMatch> ZpSlowMatchPtr;
typedef std::vector<ZpSlowMatchPtr> ZpSlowMatchPtrVec;

class ZeroPtmSlowMatch {
 public:
  ZeroPtmSlowMatch(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                   ZpFastMatchPtr fast_match_ptr,
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
