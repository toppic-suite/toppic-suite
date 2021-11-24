//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_SEARCH_ONE_PTM_SEARCH_PS_ALIGN_HPP_
#define TOPPIC_SEARCH_ONE_PTM_SEARCH_PS_ALIGN_HPP_

#include <vector>

#include "prsm/prsm.hpp"
#include "search/diag/diag_header.hpp"
#include "search/diag/diagonal.hpp"
#include "search/oneptmsearch/dp_pair.hpp"
#include "search/oneptmsearch/ps_align_para.hpp"

namespace toppic {

class PsAlign {
 public:
  PsAlign(const std::vector<double> &ms_masses,
          const std::vector<double> &seq_masses,
          const DiagonalPtrVec &diagonal_ptrs,
          PsAlignParaPtr para_ptr);

  void compute(ProteoformTypePtr type_ptr);

  void initDpPair();

  void dp(ProteoformTypePtr type_ptr);

  void backtrace();

  DiagHeaderPtrVec backtrace(int s);

  double getAlignScr(int s) {return align_scores_[s];}

  DiagHeaderPtrVec getDiagHeaders(int s) {return backtrack_diagonal_ptrs_[s];}

  PrsmPtr geneResult(int shift_num, ProteoformPtr proteo_ptr, DeconvMsPtrVec &deconv_ms_ptr_vec,
                     ExtendMsPtrVec &ms_three_ptr_vec, PrsmParaPtr prsm_para_ptr);

 protected:
  PsAlignParaPtr para_ptr_;;

  std::vector<double> ms_masses_;

  std::vector<double> seq_masses_;

  DiagonalPtrVec diagonal_ptrs_;

  std::vector<std::vector<int>> idxes_;

  std::vector<std::vector<bool>> penalties_;

  std::vector<DpPairPtrVec> dp_2d_pair_ptrs_;

  DpPairPtr first_pair_ptr_;

  DpPairPtr last_pair_ptr_;

  DpPairPtrVec segment_bgn_pair_ptrs_;

  DpPairPtrVec segment_end_pair_ptrs_;

  DpPairPtrVec dp_pair_ptrs_;

  DiagHeaderPtrVec2D backtrack_diagonal_ptrs_;

  std::vector<double> align_scores_;

  void dpPrep();

  DpPairPtr getTruncPre(DpPairPtr cur_pair_ptr, int s, ProteoformTypePtr type_ptr);

  DpPairPtr getShiftPre(int p, int s);
};

typedef std::shared_ptr<PsAlign> PsAlignPtr;
typedef std::vector<PsAlignPtr> PsAlignPtrVec;

} /* namespace toppic */

#endif /* PS_ALIGN_HPP_ */
