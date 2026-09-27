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

#ifndef TOPPIC_SEARCH_PTMSEARCH_MULTI_TYPE_DP_PAIR_HPP_
#define TOPPIC_SEARCH_PTMSEARCH_MULTI_TYPE_DP_PAIR_HPP_

#include <memory>
#include <vector>

#include "search/diag/diag_header.hpp"
#include "search/diag/pair.hpp"

namespace toppic {

class MultiTypeDpPair;
using MultiTypeDpPairPtr = std::shared_ptr<MultiTypeDpPair>;
using MultiTypeDpPairPtrVec = std::vector<MultiTypeDpPairPtr>;

// The number of proteoform types (COMPLETE, PREFIX, SUFFIX, INTERNAL,
// common_base/proteoform_type.hpp) a MultiTypeDpPair keeps a table row for.
constexpr int kProteoformTypeNum = 4;

// Plays the same role as DpPair (search/oneptmsearch/dp_pair.hpp) in the
// dynamic program of PsAlign, but for PsAlignMultiType
// (ps_align_multi_type.hpp): the (score, path type, previous pair) table of
// a DpPair has one row per unknown-shift count; a MultiTypeDpPair has one
// such row per (proteoform type, unknown-shift count), indexed by
// ProteoformType::getId() (0..3), so PsAlignMultiType::dp_ can fill in the
// table for every requested type in a single pass over the diagonal pairs,
// instead of PsAlign's one pass per type. Only the rows of the types a
// given alignment actually requests are ever written or read.
class MultiTypeDpPair : public Pair {
 public:
  MultiTypeDpPair(int x, int y, double pair_score, double diff, int order,
                  int n_shift, const DiagHeaderPtr& header_ptr);

  MultiTypeDpPairPtr getDiagPrevPairPtr() { return diag_prev_pair_ptr_; }

  void setDiagPrevPairPtr(const MultiTypeDpPairPtr& diag_prev_pair_ptr) {
    diag_prev_pair_ptr_ = diag_prev_pair_ptr;
  }

  double getDiff() { return diff_; }

  DiagHeaderPtr getDiagHeader() { return header_ptr_; }

  int getDiagOrder() { return order_; }

  double getPairScore() { return pair_score_; }

  double getScore(int type_id, int s) { return scores_[type_id][s]; }

  MultiTypeDpPairPtr getPrevPairPtr(int type_id, int s) {
    return prev_pair_ptrs_[type_id][s];
  }

  int getType(int type_id, int s) { return path_types_[type_id][s]; }

  void updateTable(int type_id, int s, double score, int path_type,
                   const MultiTypeDpPairPtr& prev_pair);

 private:
  DiagHeaderPtr header_ptr_;
  double diff_;
  double pair_score_;
  int order_;
  MultiTypeDpPairPtr diag_prev_pair_ptr_;
  // [type_id][s], type_id = ProteoformType::getId()
  std::vector<std::vector<MultiTypeDpPairPtr>> prev_pair_ptrs_;
  std::vector<std::vector<double>> scores_;
  std::vector<std::vector<int>> path_types_;
};

} /* namespace toppic */

#endif
