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


#ifndef PROT_SUFFIX_SUFFIX_TREE_HPP
#define PROT_SUFFIX_SUFFIX_TREE_HPP

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "protein_db.hpp"
#include "suffix_position.hpp"

namespace prot {

namespace suffix {

class Edge;
typedef std::shared_ptr<Edge> EdgePtr;

class Suffix;
typedef std::shared_ptr<Suffix> SuffixPtr;

class Node;
typedef std::shared_ptr<Node> NodePtr;

class SuffixTree;

class SuffixTree : public std::enable_shared_from_this<SuffixTree> {
 public:
  SuffixTree(std::string text, ProteinDBPtr database):
      seq_idx_(0),
      leaf_created_this_step_(0),
      text_(text),
      database_(database) {}

  void init();

  NodePtr getRoot();

  int getSeqIndex() {return seq_idx_;}

  int getLeafCreatedThisStep() {return leaf_created_this_step_;}

  void increaseLeafCreated() {leaf_created_this_step_++;}

  std::vector<SuffixPosPtr> search(std::string target);

  EdgePtr findMatchEdge(std::string target);

  char charAt(size_t idx) {return text_.at(idx);}

 private:
  NodePtr root;

  int seq_idx_;

  int leaf_created_this_step_;

  std::string text_;

  ProteinDBPtr database_;

  void addPrefix(SuffixPtr active, int endIndex);

  void updateSuffixNode(NodePtr node, NodePtr suffixNode);

  void determineSuffixPos();
};

typedef std::shared_ptr<SuffixTree> SuffixTreePtr;

}  // namespace suffix

}  // namespace prot
#endif
