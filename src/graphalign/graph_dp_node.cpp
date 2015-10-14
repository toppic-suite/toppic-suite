#include <limits>

#include "graph_align/graph_dp_node.hpp"

namespace prot {

GraphDpNode::GraphDpNode(int first_idx, int second_idx, 
                         double node_score, int n_shift) {
  first_idx_ = first_idx;
  second_idx_ = second_idx;
  node_score_ = node_score;
  for(int i=0;i<n_shift+1;i++){
    prev_node_ptrs_.push_back(nullptr);
    scores_.push_back(-std::numeric_limits<double>::max());
    types_.push_back(GRAPH_ALIGN_TYPE_NULL);
    best_node_ptrs_.push_back(nullptr);
    best_node_scores_.push_back(-std::numeric_limits<double>::max());
  }
}

void GraphDpNode::updateTable(int s,double score,int path_type, 
                              GraphDpNodePtr prev_node_ptr) {
    scores_[s] = score;
    types_[s] = path_type;
    prev_node_ptrs_[s] = prev_node_ptr;
}

void GraphDpNode::updateBestNode(int s, double score, 
                                 GraphDpNodePtr best_node_ptr) {
  best_node_scores_[s] = score;
  best_node_ptrs_[s] = best_node_ptr;
}

} /* namespace prot */
