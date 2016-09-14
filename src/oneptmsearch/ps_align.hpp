// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#ifndef PROT_PS_ALIGN_HPP_
#define PROT_PS_ALIGN_HPP_

#include "prsm/prsm.hpp"
#include "oneptmsearch/ps_align_para.hpp"
#include "oneptmsearch/dp_pair.hpp"
#include "oneptmsearch/basic_diag_pair.hpp"

namespace prot {

class PSAlign {
 public:
  PSAlign();
  PSAlign(const std::vector<double> &ms_masses,
          const std::vector<double> &seq_masses,
          const BasicDiagonalPtrVec &diagonal_ptrs,
          PsAlignParaPtr para_ptr);
  void compute(AlignTypePtr type_ptr);
  void initDPPair();
  void dp(AlignTypePtr type_ptr);
  void backtrace();
  DiagonalHeaderPtrVec backtrace(int s);

  double getAlignScr(int s){return align_scores_[s];};
  DiagonalHeaderPtrVec getDiagonalHeaders(int s){return backtrack_diagonal_ptrs_[s];};

  PrsmPtr geneResult(int shift_num, ProteoformPtr proteo_ptr, DeconvMsPtrVec &deconv_ms_ptr_vec,
          ExtendMsPtrVec &ms_three_ptr_vec, PrsmParaPtr prsm_para_ptr);

 protected:
  PsAlignParaPtr para_ptr_;;
  std::vector<double> ms_masses_;
  std::vector<double> seq_masses_;
  BasicDiagonalPtrVec diagonal_ptrs_;
  std::vector<std::vector<int>> idxes_;
  std::vector<std::vector<bool>> penalties_;

  std::vector<DPPairPtrVec> dp_2d_pair_ptrs_;
  DPPairPtr first_pair_ptr_;
  DPPairPtr last_pair_ptr_;
  DPPairPtrVec segment_bgn_pair_ptrs_;
  DPPairPtrVec segment_end_pair_ptrs_;
  DPPairPtrVec dp_pair_ptrs_;
  DiagonalHeaderPtrVec2D backtrack_diagonal_ptrs_;
  std::vector<double> align_scores_;

  void dpPrep();
  DPPairPtr getTruncPre(DPPairPtr cur_pair_ptr,int s, AlignTypePtr type_ptr);
  DPPairPtr getShiftPre(DPPairPtr cur_pair_ptr,int p,int s,AlignTypePtr type_ptr);
};

typedef std::shared_ptr<PSAlign> PSAlignPtr;
typedef std::vector<PSAlignPtr> PSAlignPtrVec;
} /* namespace prot */

#endif /* PS_ALIGN_HPP_ */
