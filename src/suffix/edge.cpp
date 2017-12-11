//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

namespace prot {
namespace suffix {

Edge::Edge(int beginIndex, int endIndex, Node *startNode) {
  this->beginIndex = beginIndex;
  this->endIndex = endIndex;
  this->startNode = startNode;
  this->endNode = new Node(startNode->getSuffixTree(), NULL);
  isleaf = false;
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
Node * Edge::splitEdge(Suffix * suffix) {
  remove();
  Node *breakNode = new Node(suffix->getOriginNode()->getSuffixTree(), NULL);
  Edge *newEdge = new Edge(beginIndex, beginIndex + suffix->getSpan(),
                           suffix->getOriginNode(), breakNode);
  newEdge->insert();
  breakNode->setSuffixNode(suffix->getOriginNode());
  beginIndex += suffix->getSpan() + 1;
  startNode = breakNode;
  insert();
  return breakNode;
}

/**
 * insert the edge to the list associated with the startNode
 * Each (internal) node maintains a set of edges starting from it
 * */
void Edge::insert() {
  startNode->addEdge(beginIndex, this);
}

/**
 * remove the edge from the set of edges associated with the startNode
 * */
void Edge::remove() {
  startNode->removeEdge(beginIndex);
}

/**
 * get the startNode (incoming node) of an edge
 * */
Node * Edge::getStartNode() {
  return startNode;
}

/**
 * set the startNode of an edge
 * */
void Edge::setStartNode(Node * startNode) {
  this->startNode = startNode;
}

/**
 * get the endNode of an edge
 * */
Node * Edge::getEndNode() {
  return endNode;
}

/**
 * newly added, test if an edge has an edge node (internal node if yes)
 * */
bool Edge::hasEndNode() {
  return endNode != NULL;
}

/**
 * newly added, get an element on the edge
 * */
char Edge::getItemAt(int j) {
  return startNode->getSuffixTree()->text.at(j);
}

void Edge::setLeafEdge() {
  isleaf = true;
  leaf = new LeafEdge(beginIndex,  endIndex,  startNode);
}
}  // namespace suffix
}  // namespace prot
