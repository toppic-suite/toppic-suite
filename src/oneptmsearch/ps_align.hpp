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


#ifndef PROT_PS_ALIGN_HPP_
#define PROT_PS_ALIGN_HPP_

#include <vector>

#include "prsm/prsm.hpp"
#include "oneptmsearch/ps_align_para.hpp"
#include "oneptmsearch/dp_pair.hpp"
#include "oneptmsearch/basic_diag_pair.hpp"

namespace prot {

class PSAlign {
 public:
  PSAlign() {}

  PSAlign(const std::vector<double> &ms_masses,
          const std::vector<double> &seq_masses,
          const BasicDiagonalPtrVec &diagonal_ptrs,
          PsAlignParaPtr para_ptr);

  void compute(AlignTypePtr type_ptr);

  void initDPPair();

  void dp(AlignTypePtr type_ptr);

  void backtrace();

  DiagonalHeaderPtrVec backtrace(int s);

  double getAlignScr(int s) {return align_scores_[s];}

  DiagonalHeaderPtrVec getDiagonalHeaders(int s) {return backtrack_diagonal_ptrs_[s];}

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

  DPPairPtr getTruncPre(DPPairPtr cur_pair_ptr, int s, AlignTypePtr type_ptr);
  DPPairPtr getShiftPre(DPPairPtr cur_pair_ptr, int p,int s,AlignTypePtr type_ptr);
};

typedef std::shared_ptr<PSAlign> PSAlignPtr;
typedef std::vector<PSAlignPtr> PSAlignPtrVec;
} /* namespace prot */

#endif /* PS_ALIGN_HPP_ */
