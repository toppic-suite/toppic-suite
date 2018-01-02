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


#ifndef PROT_SUFFIX_NODE_HPP
#define PROT_SUFFIX_NODE_HPP

#include <vector>

#include "edge.hpp"
#include "suffix_tree.hpp"

namespace prot {

namespace suffix {

class Node;
typedef std::shared_ptr<Node> NodePtr;

class Node {
 public:
  Node(SuffixTreePtr suffixTree, NodePtr suffixNode);

  char charAt(int index);

  void addEdge(int charIndex, EdgePtr edge) {
    edges_[getCharacterIndex(charAt(charIndex))] = edge;
  }

  void removeEdge(int charIndex) {
    edges_[getCharacterIndex(charAt(charIndex))] = nullptr;
  }

  void removeEdge(char ch) {
    edges_[getCharacterIndex(ch)] = nullptr;
  }

  EdgePtr findEdge(char ch) {return edges_[getCharacterIndex(ch)];}

  EdgePtr getEdge(int index);

  int getCharacterIndex(char ch);

  NodePtr getSuffixNode() {return suffixNode;}

  void setSuffixNode(NodePtr suffixNode) {this->suffixNode = suffixNode;}

  bool hasSuffixNode() {return suffixNode != nullptr;}

  SuffixTreePtr getSuffixTree() {return suffixTree;}

 private:
  SuffixTreePtr suffixTree;

  NodePtr suffixNode;

  std::vector<EdgePtr> edges_;
};

}  // namespace suffix

}  // namespace prot
#endif
