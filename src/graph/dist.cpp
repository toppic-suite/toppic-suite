// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include "dist.hpp"

namespace prot {

int getVecIndex(int v1, int v2, int gap) {
  int index = gap * v1 + (v2 - v1);
  return index;
}

void addToDistVec(MassGraphPtr graph_ptr, const std::vector<std::vector<std::set<int>>> & dist_vecs,
                  int node_num, int mod_num, DistVec & dist_vec, int gap) {
  std::set<Dist> dist_set;
  for (int i = 0; i < node_num - 1; i++) {
    for (int j = i + 1; j < node_num && j <= i + gap; j++) {
      int index = getVecIndex(i, j, gap);
      for (std::set<int>::iterator it=dist_vecs[index][mod_num].begin();
           it != dist_vecs[index][mod_num].end(); it++) {
        if (*it == 0) continue;
        Dist tmp = Dist(graph_ptr, *it, i, j);
        auto search = dist_set.find(tmp);
        if (search != dist_set.end()) {
          search->pair_ij_.push_back(std::pair<int, int>(i, j));
        } else {
          dist_set.insert(tmp);
        }
      }
    }
  }

  std::copy(dist_set.begin(), dist_set.end(), std::back_inserter(dist_vec));
}

}  // namespace prot
