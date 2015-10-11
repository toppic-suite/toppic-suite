#ifndef PROT_GRAPH_RESULT_NODE_HPP_
#define PROT_GRAPH_RESULT_NODE_HPP_

#include <memory>
#include <vector>

namespace prot {

class GraphResultNode { 
 public:
  GraphResultNode(GraphDpNodePtr node_ptr, int shift_num) {
    first_idx_ = node_ptr->getFirstIdx();
    second_idx_ = node_ptr->getSecondIdx();
    type_ = node_ptr->getType(shift_num);
  }

  int getFirstIdx() {return first_idx_;}
  int getSecondIdx() {return second_idx_;}
  int getType(){return type_;}

 private:
  int first_idx_;
  int second_idx_;
  int type_;
};

typedef std::shared_ptr<GraphResultNode> GraphResultNodePtr;
typedef std::vector<GraphResultNodePtr> GraphResultNodePtrVec;
typedef std::vector<GraphResultNodePtrVec> GraphResultNodePtrVec2D;

} /* namespace prot */

#endif /* GRAPH_RESULT_NODE_HPP_ */
