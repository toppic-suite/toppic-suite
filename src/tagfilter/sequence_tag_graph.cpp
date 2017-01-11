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

#include <algorithm>

#include "sequence_tag_graph.hpp"

namespace prot{

void SequenceTagGraph::SetAminoAcidsArray(std::vector<std::pair<double, std::string> > aa) {
  AminoAcidsArray = aa;
  MaxAminoAcidMass = -1;
  MinAminoAcidMass = 1e6;
  for (size_t i = 0; i < aa.size(); i++) {
    if (aa[i].first > MaxAminoAcidMass) MaxAminoAcidMass = aa[i].first;
    if (aa[i].first < MinAminoAcidMass) MinAminoAcidMass = aa[i].first;
  }
  //std::cout << "MaxAminoAcidMass " << MaxAminoAcidMass << " MinAminoAcidMass " << MinAminoAcidMass << std::endl;
}

void SequenceTagGraph::CollectSequenceTagGraphEdges() {
  for (size_t i = 0; i < PeakList.size(); i++) {
    for (size_t j = i + 1; j < PeakList.size(); j++) {
      double mass_gap = PeakList[j]->getMass() - PeakList[i]->getMass();
      double mass_err = (PeakList[j]->getMass() + PeakList[i]->getMass()) * err_tolerance / 2;
      double max_mass_gap = mass_gap + mass_err;
      double min_mass_gap = mass_gap - mass_err;

      GraphEdgePtr peakGap = std::make_shared<GraphEdge>(i, j, mass_gap);
      if (min_mass_gap > MaxAminoAcidMass) break;
      if (max_mass_gap < MinAminoAcidMass) continue;
      for (std::pair<double, std::string> aa : AminoAcidsArray) {
        if (min_mass_gap < aa.first && max_mass_gap > aa.first) {
          peakGap->AddMatchedAminoAcid(aa);
        }
      }
      if (peakGap->AminoAcidList.size() > 0) {
        AddEdge(peakGap);
      }
    }
  }
}

std::vector<std::string> SequenceTagGraph::FindSequenceTags() {
  std::vector<std::vector<int> > componentSet = ConnnectedComponents();
  std::vector<int> startNodeSet;
  for (std::vector<int> comp : componentSet) {
    if (static_cast<int>(comp.size()) < MinTaglen) continue;
    startNodeSet.insert(startNodeSet.end(), comp.begin(), comp.end());
  }

  for (int node : startNodeSet) {
    std::vector<GraphEdgePtr> edges = OutEdges(node);
    if (edges.size() < 1) continue;

    FindPaths(node);
  }

  std::vector<std::string> seq_tag_list;
  std::copy(SeqTagSet.begin(), SeqTagSet.end(), std::back_inserter(seq_tag_list));
  std::sort(seq_tag_list.begin(), seq_tag_list.end());
  return seq_tag_list;
}

std::vector<std::vector<int> > SequenceTagGraph::ConnnectedComponents() {
  std::vector<std::vector<int> > componentSet;
  std::vector<bool> visited(GetNodeCount());
  std::fill(visited.begin(), visited.end(), false); 

  for (int i = 0; i < GetNodeCount(); i++) {
    if (visited[i]) continue;

    std::vector<int> component;
    std::deque<int> neighbors;
    neighbors.push_back(i);

    while (true) {
      if (neighbors.size() < 1) break;

      int j = neighbors.front();
      neighbors.pop_front();;

      if (visited[j]) continue;

      visited[j] = true;

      component.push_back(j);

      for (GraphEdgePtr edge : _adjList[j]) {
        if (visited[edge->Node2]) continue;
        neighbors.push_back(edge->Node2);
      }

    }
    componentSet.push_back(component);
  }

  return componentSet;
}

void SequenceTagGraph::FindPaths(int node, bool firstCall, GraphEdgePtr e) {
  if (firstCall) {
    NodeVisitFlag.resize(_adjList.size());
    std::fill(NodeVisitFlag.begin(), NodeVisitFlag.end(), false);
    EdgeList.clear();
  } else {
    if (static_cast<int>(EdgeList.size()) >= MaxTagLen) {
      // need to reverse the EdgeList
      std::vector<GraphEdgePtr> tmp = EdgeList;
      std::reverse(tmp.begin(), tmp.end());
      ProcessPath(tmp); 
      // just clear
      EdgeList.clear();
    }
  }
  if (e != nullptr) EdgeList.push_back(e);

  NodeVisitFlag[node] = true;

  bool flag = false;
  for (GraphEdgePtr edge : OutEdges(node)) {
    if (!NodeVisitFlag[edge->Node2]) {
      flag = true;
      FindPaths(edge->Node2, false, edge);
    }
  }

  if (!flag) {
    // we need to reverse EdgeList
    //    std::cout << EdgeList[0]->Node1 << " " << EdgeList[0]->Node2 << std::endl;
    std::vector<GraphEdgePtr> tmp = EdgeList;
    std::reverse(tmp.begin(), tmp.end());

    bool t = ProcessPath(tmp);
    if (t == false) return;
  }
  NodeVisitFlag[node] = false;

  if (EdgeList.size() > 0) EdgeList.pop_back();
}

bool SequenceTagGraph::ProcessPath(std::vector<GraphEdgePtr> EdgeList) {
  if (static_cast<int>(EdgeList.size()) < MinTaglen ) {
    return true;
  }

  bool added = false;
  std::vector<std::string> seq_tag_list = EnumerateAllSequenceTags(EdgeList);
  for (std::string segTag : seq_tag_list) {
    if (SeqTagSet.find(segTag) == SeqTagSet.end()) {
      SeqTagSet.insert(segTag);
      added = true;
    }
  }

  return true;
}

std::vector<std::string> SequenceTagGraph::EnumerateAllSequenceTags(std::vector<GraphEdgePtr> edges) {
  std::reverse(edges.begin(), edges.end());
  std::vector<std::string> seq_tag_list;
  int tagLength = static_cast<int>(edges.size());
  std::vector<std::vector<std::pair<double, std::string> > > listOfAminoAcids;

  for (int j = 0; j < tagLength; j++) {
    listOfAminoAcids.push_back(edges[j]->AminoAcidList);
  }

  std::vector<int> indexArray(tagLength);
  std::fill(indexArray.begin(), indexArray.end(), 0);
  size_t totalCombinations = 1;

  for (size_t i = 0; i < listOfAminoAcids.size(); i++) {
    totalCombinations *= listOfAminoAcids[i].size();
  }

  for (size_t e = 0; e < totalCombinations; e++) {
    std::string tag_str = "";
    double mass = 0.0;
    for (size_t i = 0; i < indexArray.size(); i++) {
      tag_str += listOfAminoAcids[i][indexArray[i]].second;
      mass += listOfAminoAcids[i][indexArray[i]].first; 
    }
    double mass_gap = PeakList[edges[indexArray.size() - 1]->Node2]->getMass() -
        PeakList[edges[0]->Node1]->getMass();
    double mass_err = fabs(mass_gap - mass);
    double mass_th = mass_gap * err_tolerance;

    if (mass_err < mass_th) {
      //std::cout << "tag_str " << tag_str << std::endl;
      seq_tag_list.push_back(tag_str);
      std::reverse(tag_str.begin(), tag_str.end());
      //std::cout << "r_tag_str " << tag_str << std::endl;
      seq_tag_list.push_back(tag_str);
    }
    for (int i = static_cast<int>(indexArray.size()) - 1; i >= 0; i--) {
      if (indexArray[i] == static_cast<int>(listOfAminoAcids[i].size() - 1)) {
        indexArray[i] = 0;
      } else {
        indexArray[i]++;
        break;
      }
    }
  }

  return seq_tag_list;
}

}
