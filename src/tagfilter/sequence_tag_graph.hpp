// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef PROT_SEQ_TAG_GARPH
#define PROT_SEQ_TAG_GARPH

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

#include "graph_edge.hpp"
#include "peak_node.hpp"

namespace prot {

class SequenceTagGraph {
 public:
  SequenceTagGraph(int minTaglen = 5, int maxTagLen = 8):
      MinTaglen(minTaglen), MaxTagLen(maxTagLen) {}

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

}
#endif
