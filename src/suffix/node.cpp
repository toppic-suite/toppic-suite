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

namespace prot {
namespace suffix {

const int aminoAcidIndex[26] = {0, -1, 4, 3, 5, 12, 7, 8, 9, -1, 10,
  9, 11, 2, -1, 13, 6, 1, 14, 15, -1, 18, 16, -1, 17, -1};

Node::Node(SuffixTree *suffixTree, Node *suffixNode) {
  this->suffixTree = suffixTree;
  this->suffixNode = suffixNode;
  edges = new Edge*[21];
  for (int i = 0; i < 21; i++) {
    edges[i] = NULL;
  }
}

char Node::charAt(int index) {
  return suffixTree->text.at(index);
}

Node::Node(Node *node) {
  suffixNode = NULL;
  edges = NULL;
  this->suffixTree = node->suffixTree;
  this->suffixNode = NULL;
}

void Node::removeAll() {
  for (int i = 0; i < 21; i++) {
    edges[i] = NULL;
  }
  edges = NULL;
}

Edge * Node::getEdge(int index) {
  if (index < 0 || index > 20) {
    std::cerr << "invalid index" << std::endl;
    exit(0);
  }

  return edges[index];
}

// get the index corresponding to one of the 22 characters
int Node::getCharacterIndex(char ch) {
  if (ch == '#')
    return 19;
  else if (ch == '$')
    return 20;
  else
    return aminoAcidIndex[ch-'A'];
}
}  // namespace suffix
}  // namespace prot
