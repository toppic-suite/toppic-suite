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


#ifndef PROT_SUFFIX_LEAFEDGE_HPP
#define PROT_SUFFIX_LEAFEDGE_HPP

#include "edge.hpp"
#include "node.hpp"
#include "suffix_position.hpp"

namespace prot {
namespace suffix {

class Suffix;
class SuffixTree;

class LeafEdge: public Edge {
 public:
  LeafEdge(int beginIndex, int endIndex, Node * startNode):
      Edge(beginIndex, endIndex, startNode),
      leafIndex(-1),
      seqNum(-1) {}

  LeafEdge(int beginIndex, int endIndex, Node *startNode, Node *endNode):
      Edge(beginIndex, endIndex, startNode, endNode),
      leafIndex(-1),
      seqNum(-1) {}

  int getLeafIndex() {return leafIndex;}

  void setLeafIndex(int index) {leafIndex = index;}

  int getSeqNum() {return seqNum;}

  void setSeqNum(int num) {seqNum = num;}

  SuffixPosition * getSuffixPosition();

  int leafIndex;

 private:
  int seqNum;
};
}  // namespace suffix
}  // namespace prot
#endif
