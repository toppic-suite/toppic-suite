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


#ifndef PROT_SUFFIX_EDGE_HPP
#define PROT_SUFFIX_EDGE_HPP

#include <string>

namespace prot {
namespace suffix {

class Node;
class LeafEdge;
class Suffix;
class Edge {
 public:
  Edge(int beginIndex, int endIndex, Node *startNode);

  Edge(int beginIndex, int endIndex, Node *startNode, Node *endNode):
    beginIndex(beginIndex),
    endIndex(endIndex),
    startNode(startNode),
    endNode(endNode) {isleaf = false;}

  Node *splitEdge(Suffix *suffix);

  void insert();

  void remove();

  int getSpan() {return endIndex - beginIndex;}

  int getLength() {return endIndex - beginIndex + 1;}

  int getBeginIndex() {return beginIndex;}

  int getEndIndex() {return endIndex;}

  void setEndIndex(int endIndex) {this->endIndex = endIndex;}

  Node * getStartNode();

  void setStartNode(Node * startNode);

  Node * getEndNode();

  bool hasEndNode();

  char getItemAt(int j);

  bool isleaf;

  LeafEdge * leaf;

  void setLeafEdge();

 protected:
  int beginIndex;
  int endIndex;
  int span;
  Node * startNode;
  Node * endNode;
};
}  // namespace suffix
}  // namespace prot
#endif
