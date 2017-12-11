// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#ifndef PROT_TAG_VAR_SEQ_TAG_GARPH
#define PROT_TAG_VAR_SEQ_TAG_GARPH

#include <vector>
#include <queue>
#include <stack>
#include <deque>
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <deque>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <string>

#include "taglong/peak_node.hpp"
#include "graph_edge.hpp"

namespace prot {

class SequenceTagGraph {
 public:
  explicit SequenceTagGraph(int minTaglen = 5, int maxTagLen = 8):
      MinTaglen(minTaglen),
      MaxTagLen(maxTagLen) {}

  void SetNodeCount(size_t nodeCount) {
    _adjList.resize(nodeCount);
    _hasInEdge.resize(nodeCount);
  }

  void AddEdge(GraphEdgePtr edge) {
    _adjList[edge->Node1].push_back(edge);
    _hasInEdge[edge->Node2] = true;
  }

  void SetPeakList(std::vector<PeakNodePtr> p) {PeakList = p;}

  void SetAminoAcidsArray(std::vector<std::pair<double, std::string> > aa);

  bool HasInEdge(int node) { return _hasInEdge[node];}

  std::vector<GraphEdgePtr> OutEdges(int node) {return _adjList[node];}

  int GetNodeCount() { return static_cast<int>(_adjList.size());}

  std::vector<std::vector<int> > ConnnectedComponents();

  void CollectSequenceTagGraphEdges();

  void FindPaths(int node, bool firstCall = true, GraphEdgePtr e = nullptr);

  bool ProcessPath(std::vector<GraphEdgePtr> EdgeList);

  std::vector<std::string> EnumerateAllSequenceTags(std::vector<GraphEdgePtr> edges);

  std::vector<std::string> FindSequenceTags();

  void setErrTolerance(double err) {err_tolerance = err;}

 private:
  int MinTaglen, MaxTagLen;

  double MaxAminoAcidMass, MinAminoAcidMass;

  double err_tolerance;

  std::vector<std::vector<GraphEdgePtr> > _adjList;

  std::vector<bool> _hasInEdge;

  std::vector<GraphEdgePtr> EdgeList;

  std::vector<bool> NodeVisitFlag;

  std::vector<PeakNodePtr> PeakList;

  std::vector<std::pair<double, std::string> > AminoAcidsArray;

  std::unordered_set<std::string> SeqTagSet;
};

}  // namespace prot

#endif
