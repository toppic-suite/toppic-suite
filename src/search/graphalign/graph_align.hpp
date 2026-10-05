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

#ifndef TOPPIC_SEARCH_GRAPHALIGN_GRAPH_ALIGN_HPP_
#define TOPPIC_SEARCH_GRAPHALIGN_GRAPH_ALIGN_HPP_

#include "prsm/prsm.hpp"
#include "search/diag/diag_header.hpp"
#include "search/graph/graph.hpp"
#include "search/graph/proteo_graph.hpp"
#include "search/graph/spec_graph.hpp"
#include "search/graphalign/graph_align_mng.hpp"
#include "search/graphalign/graph_result_node.hpp"
#include "seq/fasta_seq.hpp"

namespace toppic {

class GraphAlign {
 public:
  GraphAlign(const GraphAlignMngPtr& mng_ptr,
             const ProteoGraphPtr& proteo_graph_ptr,
             const SpecGraphPtr& spec_graph_ptr,
             const FastaSeqPtr& fasta_seq_ptr);

  void process();

  PrsmPtr geneResult(int s, int m);

  PrsmPtr geneResult(int s);

 private:
  GraphAlignMngPtr mng_ptr_;

  ProteoGraphPtr proteo_graph_ptr_;

  MassGraphPtr pg_;

  FastaSeqPtr fasta_seq_ptr_;

  int proteo_ver_num_;

  SpecGraphPtr spec_graph_ptr_;

  MassGraphPtr sg_;

  int spec_ver_num_;

  // the spectrum distances within the range of the proteoform distances, as
  // indexes of the spectrum graph's distance list, in increasing order
  std::vector<int> spec_dist_idxes_;

  // for each proteoform N-terminal PTM mass of the manager, the proteoform
  // vertices an alignment may start at with it: a truncated start (protein
  // residue 3 or later) on an amino acid the modification allows
  std::vector<std::vector<char>> n_term_starts_;

  // The DP table: a cell per pair of a proteoform graph vertex i and a
  // spectrum graph vertex j, cell index i * spec_ver_num_ + j, and in each
  // cell a layer per unexpected shift number s and variable modification
  // number m, entry index (cell * shift_layer_num_ + s) * mod_layer_num_ + m.
  // A previous node is a cell index, -1 for none.
  int shift_layer_num_;
  int mod_layer_num_;
  std::vector<signed char> prev_edge_types_;
  std::vector<signed char> prev_edge_mod_nums_;
  std::vector<int> prev_nodes_;
  std::vector<double> best_scores_;
  std::vector<int> best_shift_nodes_;
  std::vector<double> best_shift_scores_;

  // The consistent pairs of a cell and a variable modification number m:
  // the cells cons_prev_cells_[cons_begins_[k]] ... [cons_begins_[k + 1] - 1]
  // with k = cell * mod_layer_num_ + m, in the order they were found.
  std::vector<int> cons_begins_;
  std::vector<int> cons_prev_cells_;

  GraphResultNodePtrVec3D result_nodes_;

  GraphResultNodePtrVec2D nodes_2d_;

  DiagHeaderPtrVec diag_headers_;

  DiagHeaderPtrVec2D diag_headers_2d_;

  size_t getCellIdx(int i, int j) {
    return static_cast<size_t>(i) * spec_ver_num_ + j;
  }

  size_t getEntryIdx(size_t cell, int s, int m) {
    return (cell * shift_layer_num_ + s) * mod_layer_num_ + m;
  }

  void getConsistentPairs();

  void initTable();

  int compBestVariableNode(size_t cell, int s, int m, int& best_edge_mod_num);

  int compBestShiftNode(int i, int j, int s, int m);

  void updateBestShiftNode(int i, int j, int s, int m);

  void dp();

  GraphResultNodePtrVec backtrace(int s, int m);

  void backtrace();

  void getNodeDiagonals(int s, int m);

  void geneHeaders();
};

using GraphAlignPtr = std::shared_ptr<GraphAlign>;

}  // namespace toppic

#endif
