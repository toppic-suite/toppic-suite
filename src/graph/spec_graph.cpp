//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#include <utility>
#include <algorithm>
#include <vector>
#include <set>

#include "common/util/logger.hpp"
#include "graph/spec_graph.hpp"

namespace toppic {

SpecGraph::SpecGraph(SpectrumSetPtr spec_set_ptr, PrmPeakPtrVec peak_vec,
                     MassGraphPtr graph_ptr, double convert_ratio) {
  spec_set_ptr_ = spec_set_ptr;
  peak_vec_ = peak_vec;
  graph_ptr_ = graph_ptr;
  node_num_ = num_vertices(*graph_ptr.get());
  pair_num_ = node_num_ * (node_num_ + 1) /2;
  compSpecDistances(convert_ratio);
}

int SpecGraph::getVecIndex(int v1, int v2) {
  int index =  (node_num_ + node_num_ - v1 + 1) * v1 /2  + (v2 - v1);
  return index;
}

int SpecGraph::getPeakDist(int v1, int v2) {
  int index = getVecIndex(v1, v2);
  return peak_dists_[index];
}

void SpecGraph::compSpecDistances(double convert_ratio) {
  std::set<Dist> dist_set;
  int count = 0;
  peak_dists_ = std::vector<int>(pair_num_, 0);
  for (size_t i = 0; i < peak_vec_.size() - 1; i++) {
    for (size_t j = i+1; j < peak_vec_.size(); j++) {
      double dist = peak_vec_[j]->getPosition() - peak_vec_[i]->getPosition();

      int int_dist = std::round(dist * convert_ratio);

      Dist tmp = Dist(graph_ptr_, int_dist, i, j);
      auto search = dist_set.find(tmp);

      if (search != dist_set.end()) {
        search->pair_ij_.push_back(std::pair<int, int>(i, j));
      } else {
        dist_set.insert(tmp);
      }

      int index = getVecIndex(i, j);
      peak_dists_[index] = int_dist;
      count++;
    }
  }
  LOG_DEBUG("count " << count);
  std::copy(dist_set.begin(), dist_set.end(), std::back_inserter(dist_));
}

}  // namespace toppic
