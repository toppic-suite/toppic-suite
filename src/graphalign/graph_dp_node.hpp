#ifndef PROT_GRAPH_DP_NODE_HPP_
#define PROT_GRAPH_DP_NODE_HPP_

#include <memory>
#include <vector>

namespace prot {

#define GRAPH_ALIGN_TYPE_NULL -1
#define GRAPH_ALIGN_TYPE_VARIABLE 0 
#define GRAPH_ALIGN_TYPE_UNEXPECTED 1

class GraphDpNode;
typedef std::shared_ptr<GraphDpNode>  GraphDpNodePtr;
typedef std::vector<GraphDpNodePtr> GraphDpNodePtrVec;
typedef std::vector<GraphDpNodePtrVec> GraphDpNodePtrVec2D;

typedef std::weak_ptr<GraphDpNode>  GraphDpNodeWeakPtr;
typedef std::vector<GraphDpNodeWeakPtr> GraphDpNodeWeakPtrVec;
typedef std::vector<GraphDpNodeWeakPtrVec> GraphDpNodeWeakPtrVec2D;

class GraphDpNode { 
 public:
  GraphDpNode(int first_idx, int second_idx, double node_score,
              int n_unknown_shifts, int max_known_mods);
  int getFirstIdx() {return first_idx_;}
  int getSecondIdx() {return second_idx_;}

  double getNodeScore() {return node_score_;}
  int getPrevEdgeType(int s, int m){return prev_edge_types_[s][m];}
  int getPrevEdgeModNum(int s, int m){return prev_edge_mod_nums_[s][m];}

  double getBestScore(int s, int m) {return best_scores_[s][m];}

  GraphDpNodePtr getPrevNodePtr(int s, int m){return prev_node_ptrs_[s][m].lock();}

  void updateTable(int s, int m, int path_type, int mod_num,
                   GraphDpNodePtr prev_node_ptr, int score);

  void updateBestShiftNode(int s, int m, double score, GraphDpNodePtr prev_node_ptr);

  double getBestShiftScore(int s, int m) {return best_shift_scores_[s][m];}

  GraphDpNodePtr getBestShiftNodePtr(int s, int m) {return best_shift_node_ptrs_[s][m].lock();}


 private:
  int first_idx_;
  int second_idx_;
  // the score for current node
  double node_score_;

  // prev edge and node for dp
  std::vector<std::vector<int>> prev_edge_types_;
  std::vector<std::vector<int>> prev_edge_mod_nums_;

  GraphDpNodeWeakPtrVec2D prev_node_ptrs_;
  std::vector<std::vector<double>> best_scores_;

  // the vector for finding shift nodes
  GraphDpNodeWeakPtrVec2D best_shift_node_ptrs_;
  std::vector<std::vector<double>> best_shift_scores_;
};


} /* namespace prot */

#endif /* GRAPH_DP_NODE_HPP_ */
