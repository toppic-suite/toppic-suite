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
    GraphDpNodePtrVec node_ptrs (max_known_mods + 1, nullptr);
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
