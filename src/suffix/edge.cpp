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


#include "node.hpp"
#include "edge.hpp"
#include "suffix.hpp"
#include "suffix_tree.hpp"

namespace toppic {
namespace suffix {

Edge::Edge(int beginIndex, int endIndex, NodePtr startNode) {
  bgn_idx_ = beginIndex;
  end_idx_ = endIndex;
  start_node_ = startNode;
  end_node_ = std::make_shared<Node>(startNode->getSuffixTree(), nullptr);
  is_leaf_ = false;
}

/**
 * When a suffix ends on an implicit node, adding a new character means I have to split an existing edge.
 * This function is called to split an edge at the point defined by the Suffix argument.
 * The existing edge loses its parent, as well as some of its leading characters.
 * The newly created edge descends from the original parent, and now has the existing edge as a child.
 *
 * Since the existing edge is getting a new parent and starting character,
 * its hash table entry will no longer be valid.  That's why it gets removed at the start of the function.
 * After the parent and start char have been recalculated, it is re-inserted.
 *
 * The number of characters stolen from the original node and given to the new node is equal to the number
 * of characters in the suffix argument, which is last - first + 1;
 */
NodePtr Edge::splitEdge(SuffixPtr suffix) {
  remove();
  NodePtr breakNode = std::make_shared<Node>(suffix->getOriginNode()->getSuffixTree(), nullptr);
  EdgePtr newEdge = std::make_shared<Edge>(bgn_idx_, bgn_idx_ + suffix->getSpan(), suffix->getOriginNode(), breakNode);
  newEdge->insert();
  breakNode->setSuffixNode(suffix->getOriginNode());
  bgn_idx_ += suffix->getSpan() + 1;
  start_node_ = breakNode;
  insert();
  return breakNode;
}

/**
 * insert the edge to the list associated with the startNode
 * Each (internal) node maintains a set of edges starting from it
 * */
void Edge::insert() {
  start_node_->addEdge(bgn_idx_, shared_from_this());
}

// remove the edge from the set of edges associated with the startNode
void Edge::remove() {
  start_node_->removeEdge(bgn_idx_);
}

// get the startNode (incoming node) of an edge
NodePtr Edge::getStartNode() {
  return start_node_;
}

// get the endNode of an edge
NodePtr Edge::getEndNode() {
  return end_node_;
}

// newly added, test if an edge has an edge node (internal node if yes)
bool Edge::hasEndNode() {
  return end_node_ != nullptr;
}

// newly added, get an element on the edge
char Edge::getItemAt(int j) {
  return start_node_->getSuffixTree()->charAt(j);
}

void Edge::setLeafEdge() {
  is_leaf_ = true;
  leaf_ = std::make_shared<LeafEdge>(bgn_idx_,  end_idx_,  start_node_);
}

}  // namespace suffix

}  // namespace toppic
