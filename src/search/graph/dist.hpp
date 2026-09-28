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

#include <utility>
#include <vector>

namespace toppic {

// The distances between the vertex pairs of a graph, one entry per distance
// in increasing order: the d-th distance is dists_[d], and its vertex pairs
// (i, j) are pairs_[begins_[d]] ... pairs_[begins_[d + 1] - 1], in the order
// they were found.
struct DistList {
  std::vector<int> dists_;
  std::vector<int> begins_{0};
  std::vector<std::pair<int, int>> pairs_;

  size_t size() const { return dists_.size(); }
};

using DistList2D = std::vector<DistList>;

// a distance between the vertices i and j of a graph
struct DistEntry {
  int dist_;
  int i_;
  int j_;
};

// Adds the entries to an empty list: one distance per value, in increasing
// order, each with its (i, j) pairs in the order of the entries. The entries
// are reordered.
void addToDistList(std::vector<DistEntry>& entries, DistList& dist_list);

// The distances of the pairs (i, j) of a proteoform graph with j - i at most
// gap and mod_num variable modifications: dist_vecs[(gap + 1) * i + j - i]
// [mod_num] holds the distances of the pair, sorted and unique, and a
// distance of 0 is left out.
void addToDistList(const std::vector<std::vector<std::vector<int>>>& dist_vecs,
                   int node_num, int mod_num, DistList& dist_list, int gap);

}  // namespace toppic

#endif
