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


#ifndef PROT_SUFFIX_TREE_HPP
#define PROT_SUFFIX_TREE_HPP

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "protein_db.hpp"
#include "suffix_position.hpp"
#include "node.hpp"

namespace prot {
namespace suffix {

class Edge;
class LeafEdge;
class Suffix;
class SuffixTree {
 public:
  explicit SuffixTree(std::string text);

  SuffixTree(std::string text, ProteinDatabase *database);

  Node * getRoot();

  int getSeqIndex() {return seqIndex;}

  int getLeafCreatedThisStep() {return leafCreatedThisStep;}

  void increaseLeafCreated() {leafCreatedThisStep++;}

  std::vector<SuffixPosition *> search(std::string target);

  Edge * findMatchEdge(std::string target);

  std::string getstring(char c);

  std::string text;

 private:
  Node * root;
  int seqIndex;
  int leafCreatedThisStep;
  ProteinDatabase * database;
  void addPrefix(Suffix * active, int endIndex);
  void updateSuffixNode(Node * node, Node * suffixNode);
  void determineSuffixPos();
};

}  // namespace suffix
}  // namespace prot
#endif
