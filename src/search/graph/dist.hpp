//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#ifndef TOPPIC_SEARCH_GRAPH_DIST_HPP_
#define TOPPIC_SEARCH_GRAPH_DIST_HPP_

#include <utility>
#include <set>
#include <vector>

#include "search/graph/graph.hpp"

namespace toppic {

class Dist{
 public:
  Dist(MassGraphPtr graph_ptr, int d, int i, int j) {
    graph_ptr_ = graph_ptr;
    dist_ = d;
    pair_ij_.push_back(std::pair<int, int>(i, j));
  }

  bool operator< (const Dist& d) const {
    return dist_ < d.dist_;
  }

  mutable std::vector<std::pair<int, int>> pair_ij_;

  int dist_;

  MassGraphPtr getGraphPtr() {return graph_ptr_;}

 private:
  MassGraphPtr graph_ptr_;
};

typedef std::vector<Dist> DistVec;
typedef std::vector<DistVec> DistVec2D;

inline bool distVecUp(const Dist & a, const Dist & b) {
  return a.dist_ < b.dist_;
}

void addToDistVec(MassGraphPtr graph_ptr,
                  const std::vector<std::vector<std::set<int>>> & dist_vecs,
                  int node_num, int mod_num, DistVec & dist_vec,
                  int gap);

}  // namespace toppic

#endif
