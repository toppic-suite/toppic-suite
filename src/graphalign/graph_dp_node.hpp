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

class GraphDpNode { 
 public:
  GraphDpNode(int first_idx, int second_idx, double node_score,
              int n_shift);
  int getFirstIdx() {return first_idx_;}
  int getSecondIdx() {return second_idx_;}

  double getScore(int s) {return scores_[s];}

  GraphDpNodePtr getPrevNodePtr(int s){return prev_node_ptrs_[s];}

  int getType(int s){return types_[s];}

  void updateTable(int s,double score,int path_type, 
                   GraphDpNodePtr prev_node_ptr);

  void updateBestNode(int s, double score, GraphDpNodePtr prev_node_ptr);

  double getBestNodeScore(int s) {return best_node_scores_[s];}

  GraphDpNodePtr getBestNodePtr(int s) {return best_node_ptrs_[s];}

  double getNodeScore() {return node_score_;}

 private:
  int first_idx_;
  int second_idx_;
  double node_score_;
  GraphDpNodePtrVec prev_node_ptrs_;
  std::vector<double> scores_;
  std::vector<int> types_;
  GraphDpNodePtrVec best_node_ptrs_;
  std::vector<double> best_node_scores_;
};

typedef std::vector<GraphDpNodePtrVec> GraphDpNodePtrVec2D;

} /* namespace prot */

#endif /* GRAPH_DP_NODE_HPP_ */
