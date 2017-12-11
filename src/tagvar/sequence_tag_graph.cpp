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

#include <algorithm>
#include <utility>
#include <deque>
#include <string>
#include <vector>

#include "base/logger.hpp"

#include "sequence_tag_graph.hpp"

namespace prot {

void SequenceTagGraph::SetAminoAcidsArray(std::vector<std::pair<double, std::string> > aa) {
  AminoAcidsArray = aa;
  MaxAminoAcidMass = -1;
  MinAminoAcidMass = 1e6;
  for (size_t i = 0; i < aa.size(); i++) {
    if (aa[i].first > MaxAminoAcidMass) MaxAminoAcidMass = aa[i].first;
    if (aa[i].first < MinAminoAcidMass) MinAminoAcidMass = aa[i].first;
  }

  LOG_DEBUG("MaxAminoAcidMass " << MaxAminoAcidMass << " MinAminoAcidMass " << MinAminoAcidMass);
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
  if (static_cast<int>(EdgeList.size()) < MinTaglen) {
    return true;
  }

  // bool added = false;
  std::vector<std::string> seq_tag_list = EnumerateAllSequenceTags(EdgeList);
  for (std::string segTag : seq_tag_list) {
    if (SeqTagSet.find(segTag) == SeqTagSet.end()) {
      SeqTagSet.insert(segTag);
      // added = true;
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
      LOG_DEBUG("tag_str " << tag_str);
      seq_tag_list.push_back(tag_str);
      std::reverse(tag_str.begin(), tag_str.end());
      LOG_DEBUG("r_tag_str " << tag_str);
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

}  // namespace prot
