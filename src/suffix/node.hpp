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


#ifndef PROT_SUFFIX_NODE_HPP
#define PROT_SUFFIX_NODE_HPP

#include "edge.hpp"
#include "leaf_edge.hpp"
#include "suffix.hpp"
#include "suffix_tree.hpp"

namespace prot {

namespace suffix {

class Node {
 public:
  Node(SuffixTree * suffixTree, Node * suffixNode);

  explicit Node(Node * node);

  char charAt(int index);

  void addEdge(int charIndex, Edge * edge) {
    edges[getCharacterIndex(charAt(charIndex))] = edge;
  }

  void removeEdge(int charIndex) {
    edges[getCharacterIndex(charAt(charIndex))] = NULL;
  }

  void removeEdge(char ch) {
    edges[getCharacterIndex(ch)] = NULL;
  }

  void removeAll();

  Edge * findEdge(char ch) {return edges[getCharacterIndex(ch)];}

  Edge * getEdge(int index);

  int getCharacterIndex(char ch);

  Node * getSuffixNode() {return suffixNode;}

  void setSuffixNode(Node * suffixNode) {this->suffixNode = suffixNode;}

  bool hasSuffixNode() {return suffixNode != NULL;}

  SuffixTree * getSuffixTree() {return suffixTree;}

 private:
  SuffixTree * suffixTree;

  Node * suffixNode;

  Edge ** edges;
};

}  // namespace suffix

}  // namespace prot
#endif
