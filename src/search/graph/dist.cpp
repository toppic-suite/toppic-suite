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

void addToDistList(std::vector<DistEntry>& entries, DistList& dist_list) {
  // a stable sort keeps the entries of one distance in their order
  std::stable_sort(entries.begin(), entries.end(),
                   [](const DistEntry& a, const DistEntry& b) {
                     return a.dist_ < b.dist_;
                   });
  dist_list.pairs_.reserve(entries.size());
  for (size_t k = 0; k < entries.size(); k++) {
    if (k == 0 || entries[k].dist_ != entries[k - 1].dist_) {
      if (k > 0) {
        dist_list.begins_.push_back(static_cast<int>(k));
      }
      dist_list.dists_.push_back(entries[k].dist_);
    }
    dist_list.pairs_.push_back(std::pair<int, int>(entries[k].i_, entries[k].j_));
  }
  if (entries.size() > 0) {
    dist_list.begins_.push_back(static_cast<int>(entries.size()));
  }
}

void addToDistList(const std::vector<std::vector<std::vector<int>>>& dist_vecs,
                   int node_num, int mod_num, DistList& dist_list, int gap) {
  std::vector<DistEntry> entries;
  for (int i = 0; i < node_num - 1; i++) {
    for (int j = i + 1; j < node_num && j <= i + gap; j++) {
      int index = getVecIndex(i, j, gap);
      const std::vector<int>& dists = dist_vecs[index][mod_num];
      for (size_t k = 0; k < dists.size(); k++) {
        if (dists[k] == 0) continue;
        entries.push_back(DistEntry{dists[k], i, j});
      }
    }
  }
  addToDistList(entries, dist_list);
}

}  // namespace toppic
