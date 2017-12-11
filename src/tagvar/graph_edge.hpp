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


#ifndef PROT_TAG_VAR_GRAPH_DEGE
#define PROT_TAG_VAR_GRAPH_DEGE

#include <vector>
#include <memory>
#include <utility>
#include <string>

namespace prot {

class GraphEdge {
 public:
  GraphEdge(int node1, int node2): Node1(node1), Node2(node2) {}

  GraphEdge(int peak1Index, int peak2Index, double peakGapDistance): Node1(peak1Index),
    Node2(peak2Index), Mass(peakGapDistance) {}

  void AddMatchedAminoAcid(std::pair<double, std::string> aa) {
    AminoAcidList.push_back(aa);
  }

  int Node1;
  int Node2;
  double Mass;

  std::vector<std::pair<double, std::string> > AminoAcidList;
};

typedef std::shared_ptr<GraphEdge> GraphEdgePtr;

}  // namespace prot

#endif
