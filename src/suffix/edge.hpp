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

#include "suffix_position.hpp"

namespace prot {

namespace suffix {

class Node;
typedef std::shared_ptr<Node> NodePtr;

class LeafEdge;
typedef std::shared_ptr<LeafEdge> LeafEdgePtr;

class Suffix;
typedef std::shared_ptr<Suffix> SuffixPtr;

class Edge;

class Edge : public std::enable_shared_from_this<Edge> {
 public:
  Edge(int beginIndex, int endIndex, NodePtr startNode);

  Edge(int beginIndex, int endIndex, NodePtr startNode, NodePtr endNode):
      bgn_idx_(beginIndex),
      end_idx_(endIndex),
      start_node_(startNode),
      end_node_(endNode),
      is_leaf_(false) {}

  NodePtr splitEdge(SuffixPtr suffix);

  void insert();

  void remove();

  int getSpan() {return end_idx_ - bgn_idx_;}

  int getBeginIndex() {return bgn_idx_;}

  int getEndIndex() {return end_idx_;}

  void setEndIndex(int endIndex) {end_idx_ = endIndex;}

  NodePtr getStartNode();

  NodePtr getEndNode();

  bool hasEndNode();

  char getItemAt(int j);

  LeafEdgePtr leaf_;

  void setLeafEdge();

 protected:
  int bgn_idx_;

  int end_idx_;

  NodePtr start_node_;

  NodePtr end_node_;

  bool is_leaf_;
};

typedef std::shared_ptr<Edge> EdgePtr;

class LeafEdge: public Edge {
 public:
  LeafEdge(int beginIndex, int endIndex, NodePtr startNode):
      Edge(beginIndex, endIndex, startNode),
      leaf_idx_(-1),
      seq_num_(-1) {}

  LeafEdge(int beginIndex, int endIndex, NodePtr startNode, NodePtr endNode):
      Edge(beginIndex, endIndex, startNode, endNode),
      leaf_idx_(-1),
      seq_num_(-1) {}

  void setLeafIndex(int index) {leaf_idx_ = index;}

  int getSeqNum() {return seq_num_;}

  void setSeqNum(int num) {seq_num_ = num;}

  SuffixPosPtr getSuffixPosition() {
    if (leaf_idx_ == -1) {
      return nullptr;
    }
    return std::make_shared<SuffixPosition>(seq_num_, leaf_idx_);
  }

 private:
  int leaf_idx_;

  int seq_num_;
};

}  // namespace suffix

}  // namespace prot

#endif
