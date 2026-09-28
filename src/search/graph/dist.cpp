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

#include "search/graph/dist.hpp"

#include <algorithm>

namespace toppic {

int getVecIndex(int v1, int v2, int gap) {
  int index = (gap + 1) * v1 + (v2 - v1);
  return index;
}

void addToDistVec(const MassGraphPtr& graph_ptr,
                  std::vector<DistEntry>& entries, DistVec& dist_vec) {
  // a stable sort keeps the entries of one distance in their order
  std::stable_sort(entries.begin(), entries.end(),
                   [](const DistEntry& a, const DistEntry& b) {
                     return a.dist_ < b.dist_;
                   });
  size_t bgn = 0;
  while (bgn < entries.size()) {
    size_t end = bgn + 1;
    while (end < entries.size() && entries[end].dist_ == entries[bgn].dist_) {
      end++;
    }
    Dist dist(graph_ptr, entries[bgn].dist_, entries[bgn].i_, entries[bgn].j_);
    dist.pair_ij_.reserve(end - bgn);
    for (size_t k = bgn + 1; k < end; k++) {
      dist.pair_ij_.push_back(std::pair<int, int>(entries[k].i_, entries[k].j_));
    }
    dist_vec.push_back(std::move(dist));
    bgn = end;
  }
}

void addToDistVec(const MassGraphPtr& graph_ptr,
                  const std::vector<std::vector<std::set<int>>>& dist_vecs,
                  int node_num, int mod_num, DistVec& dist_vec, int gap) {
  std::vector<DistEntry> entries;
  for (int i = 0; i < node_num - 1; i++) {
    for (int j = i + 1; j < node_num && j <= i + gap; j++) {
      int index = getVecIndex(i, j, gap);
      for (std::set<int>::iterator it = dist_vecs[index][mod_num].begin();
           it != dist_vecs[index][mod_num].end(); it++) {
        if (*it == 0) continue;
        entries.push_back(DistEntry{*it, i, j});
      }
    }
  }
  addToDistVec(graph_ptr, entries, dist_vec);
}

}  // namespace toppic
