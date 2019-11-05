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


#include <limits>

#include "graphalign/graph_dp_node.hpp"

namespace toppic {

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

} /* namespace toppic */
