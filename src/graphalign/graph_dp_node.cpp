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


#include <limits>

#include "graphalign/graph_dp_node.hpp"

namespace prot {

GraphDpNode::GraphDpNode(int first_idx, int second_idx, 
                         double node_score, int n_unknown_shifts, 
                         int max_known_mods) {
  first_idx_ = first_idx;
  second_idx_ = second_idx;
  node_score_ = node_score;
  for(int i=0; i < n_unknown_shifts + 1;i++){
    std::vector<int> edge_types(max_known_mods + 1, GRAPH_ALIGN_TYPE_NULL);
    prev_edge_types_.push_back(edge_types);
    std::vector<int> mod_nums (max_known_mods + 1, 0);
    prev_edge_mod_nums_.push_back(mod_nums);
    GraphDpNodePtr empty_node = nullptr;
    GraphDpNodeWeakPtrVec node_ptrs (max_known_mods + 1, empty_node);
    prev_node_ptrs_.push_back(node_ptrs);
    best_shift_node_ptrs_.push_back(node_ptrs);
    std::vector<double> scores(max_known_mods + 1, -std::numeric_limits<double>::max());
    best_scores_.push_back(scores);
    best_shift_scores_.push_back(scores);
  }
}

void GraphDpNode::updateTable(int s, int m, int path_type, int mod_num,
                              GraphDpNodePtr prev_node_ptr, int score) {
  prev_edge_types_[s][m] = path_type;
  prev_edge_mod_nums_[s][m] = mod_num;
  prev_node_ptrs_[s][m] = prev_node_ptr;
  best_scores_[s][m] = score;
}

void GraphDpNode::updateBestShiftNode(int s, int m, double score, 
                                      GraphDpNodePtr best_node_ptr) {
  best_shift_scores_[s][m] = score;
  best_shift_node_ptrs_[s][m] = best_node_ptr;
}

} /* namespace prot */
