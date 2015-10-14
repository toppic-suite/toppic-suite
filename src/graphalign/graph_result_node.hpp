#ifndef PROT_GRAPH_RESULT_NODE_HPP_
#define PROT_GRAPH_RESULT_NODE_HPP_

#include <memory>
#include <vector>

namespace prot {

class GraphResultNode { 
 public:
  GraphResultNode(GraphDpNodePtr node_ptr, int s, int m) {
    first_idx_ = node_ptr->getFirstIdx();
    second_idx_ = node_ptr->getSecondIdx();
    prev_edge_type_ = node_ptr->getPrevEdgeType(s, m);
    prev_edge_mod_num_ = node_ptr->getPrevEdgeModNum(s, m);
  }

  int getFirstIdx() {return first_idx_;}
  int getSecondIdx() {return second_idx_;}
  int getPrevEdgeType(){return prev_edge_type_;}
  int getPrevEdgeModNum(){return prev_edge_mod_num_;}

 private:
  int first_idx_;
  int second_idx_;
  int prev_edge_type_;
  int prev_edge_mod_num_;
};

typedef std::shared_ptr<GraphResultNode> GraphResultNodePtr;
typedef std::vector<GraphResultNodePtr> GraphResultNodePtrVec;
typedef std::vector<GraphResultNodePtrVec> GraphResultNodePtrVec2D;

} /* namespace prot */

#endif /* GRAPH_RESULT_NODE_HPP_ */
