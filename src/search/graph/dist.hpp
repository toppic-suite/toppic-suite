// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_SEARCH_GRAPH_DIST_HPP_
#define TOPPIC_SEARCH_GRAPH_DIST_HPP_

#include <set>
#include <utility>
#include <vector>

#include "search/graph/graph.hpp"

namespace toppic {

class Dist {
 public:
  Dist(const MassGraphPtr& graph_ptr, int d, int i, int j) {
    graph_ptr_ = graph_ptr;
    dist_ = d;
    pair_ij_.push_back(std::pair<int, int>(i, j));
  }

  bool operator<(const Dist& d) const { return dist_ < d.dist_; }

  mutable std::vector<std::pair<int, int>> pair_ij_;

  int dist_;

  MassGraphPtr getGraphPtr() { return graph_ptr_; }

 private:
  MassGraphPtr graph_ptr_;
};

using DistVec = std::vector<Dist>;
using DistVec2D = std::vector<DistVec>;

// a distance between the vertices i and j of a graph
struct DistEntry {
  int dist_;
  int i_;
  int j_;
};

// Appends to dist_vec one Dist per distance of the entries, in increasing
// order of distance, each holding its (i, j) pairs in the order of the
// entries: what a std::set<Dist> ordered by distance gives when the entries
// are added to it one by one, the pair of an entry whose distance is already
// there appended to that Dist. The entries are reordered.
void addToDistVec(const MassGraphPtr& graph_ptr,
                  std::vector<DistEntry>& entries, DistVec& dist_vec);

inline bool distVecUp(const Dist& a, const Dist& b) {
  return a.dist_ < b.dist_;
}

void addToDistVec(const MassGraphPtr& graph_ptr,
                  const std::vector<std::vector<std::set<int>>>& dist_vecs,
                  int node_num, int mod_num, DistVec& dist_vec, int gap);

}  // namespace toppic

#endif
