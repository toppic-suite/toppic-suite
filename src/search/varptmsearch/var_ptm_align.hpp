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

#ifndef TOPPIC_SEARCH_VAR_PTM_SEARCH_VAR_PTM_ALIGN_HPP_
#define TOPPIC_SEARCH_VAR_PTM_SEARCH_VAR_PTM_ALIGN_HPP_

#include <vector>

#include "prsm/prsm.hpp"
#include "search/diag/diagonal.hpp"
#include "search/varptmsearch/var_ptm_search_mng.hpp"

namespace toppic {

class VarPtmAlign {
 public:
  VarPtmAlign(const std::vector<double> &ms_masses,
              const std::vector<double> &seq_masses,
              ResSeqPtr res_seq_ptr,
              const DiagonalPtrVec &diagonal_ptrs,
              VarPtmSearchMngPtr mng_ptr);

  void initMatchTable();

  void initAllowShiftTable();

  void initScoreTable();

  void compute();

  void dp();

  DiagHeaderPtrVec backtrack(); 

  // backtrack for diagonal d and variable ptm number s
  DiagHeaderPtrVec backtrack(int d, int s); 

  /*
  double getAlignScr(int s) {return align_scores_[s];}

  DiagHeaderPtrVec getDiagHeaders(int s) {return backtrack_diagonal_ptrs_[s];}

  PrsmPtr geneResult(int shift_num, ProteoformPtr proteo_ptr, DeconvMsPtrVec &deconv_ms_ptr_vec,
                     ExtendMsPtrVec &ms_three_ptr_vec, PrsmParaPtr prsm_para_ptr);
                     */

 protected:
  std::vector<double> ms_masses_;

  std::vector<double> seq_masses_;

  ResSeqPtr res_seq_ptr_;

  DiagonalPtrVec diagonal_ptrs_;

  VarPtmSearchMngPtr mng_ptr_;

  std::vector<std::vector<int>> matches_2d_;

  std::vector<std::vector<bool>> allow_shift_2d_; 

  std::vector<std::vector<std::vector<int>>> scores_3d_;

  std::vector<std::vector<std::vector<int>>> prevs_3d_;

  /*
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
  */
};

typedef std::shared_ptr<VarPtmAlign> VarPtmAlignPtr;

} /* namespace toppic */

#endif /* PS_ALIGN_HPP_ */
