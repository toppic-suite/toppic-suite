// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_SEARCH_VARPTMSEARCH_VAR_PTM_ALIGN_HPP_
#define TOPPIC_SEARCH_VARPTMSEARCH_VAR_PTM_ALIGN_HPP_

#include <vector>

#include "prsm/prsm.hpp"
#include "search/diag/diagonal.hpp"
#include "search/varptmsearch/var_ptm_search_mng.hpp"

namespace toppic {

class VarPtmAlign {
 public:
  VarPtmAlign(const std::vector<double>& ms_masses,
              const std::vector<double>& seq_masses,
              const ResSeqPtr& sub_res_seq_ptr,
              const DiagonalPtrVec& diagonal_ptrs,
              const VarPtmSearchMngPtr& mng_ptr, int max_ptm_num);

  void initMatchTable();

  void initAllowShiftTable();

  void initScoreTable();

  void compute();

  void dp();

  void backtrack();

  // backtrack for diagonal d and variable ptm number ptm
  void backtrack(int d, int ptm);

  double getAlignScr() { return best_score_; }

  DiagHeaderPtrVec getDiagHeaders() { return backtrack_diag_header_ptrs_; }

  MassShiftPtrVec geneShiftVec(const ProteoformPtr& sub_proteo_ptr,
                               DiagHeaderPtrVec& ori_header_ptrs,
                               DiagHeaderPtrVec& refined_header_ptrs);

  // A shift whose PTMs have as many candidate residues as alterations (the
  // residues of its range that a modification of the PTM modifies, none a
  // candidate of two PTMs) has known sites: it becomes one shift per PTM on
  // its own residue, e.g. (MAQM)[Oxidation;Oxidation] -> (M)[Oxidation]AQ
  // (M)[Oxidation]. Any other shift is returned as it is.
  MassShiftPtrVec splitDeterminedShift(const MassShiftPtr& shift_ptr);

  PrsmPtr geneResult(const ProteoformPtr& sub_proteo_ptr,
                     DeconvMsPtrVec& deconv_ms_ptr_vec,
                     const PrsmParaPtr& prsm_para_ptr);

 protected:
  std::vector<double> ms_masses_;

  std::vector<double> seq_masses_;

  ResSeqPtr sub_res_seq_ptr_;

  DiagonalPtrVec diagonal_ptrs_;

  VarPtmSearchMngPtr mng_ptr_;

  std::vector<std::vector<int>> matches_2d_;

  std::vector<std::vector<bool>> allow_shift_2d_;

  std::vector<std::vector<std::vector<int>>> scores_3d_;

  std::vector<std::vector<std::vector<int>>> prevs_3d_;

  int best_d_;
  int best_ptm_num_;

  // the largest number of variable PTMs of a result: the manager's
  // var_ptm_num_, less one for a proteoform with a truncation N-terminal
  // modification
  int max_ptm_num_;
  int best_score_;

  DiagHeaderPtrVec backtrack_diag_header_ptrs_;

  ModPtrVec backtrack_mod_ptrs_;
};

using VarPtmAlignPtr = std::shared_ptr<VarPtmAlign>;

} /* namespace toppic */

#endif /* PS_ALIGN_HPP_ */
