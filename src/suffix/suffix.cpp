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
#include "leaf_edge.hpp"
#include "suffix.hpp"
#include "suffix_tree.hpp"

namespace prot {
namespace suffix {

Suffix::Suffix(Node *originNode, int beginIndex, int endIndex) {
  this->originNode = originNode;
  this->beginIndex = beginIndex;
  this->endIndex = endIndex;
}

/**
 * match the suffix along the tree from Suffix.originNode
 * */
void Suffix::canonize() {
  if (!isExplicit()) {
    Edge *edge = originNode->findEdge(originNode->charAt(beginIndex));

    int edgeSpan = edge->getSpan();
    while (edgeSpan <= getSpan()) {
      beginIndex += edgeSpan + 1;
      originNode = edge->getEndNode();
      if (beginIndex <= endIndex) {
        edge = edge->getEndNode()->findEdge(originNode->charAt(beginIndex));
        edgeSpan = edge->getSpan();
      }
    }
  }
}

/**
 * change the origin node of the suffix
 * */
void Suffix::changeOriginNode() {
  originNode = originNode->getSuffixNode();
}
}  // namespace suffix
}  // namespace prot

