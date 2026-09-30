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

#ifndef TOPPIC_TOPLIB_LIB_CLUSTER_HPP_
#define TOPPIC_TOPLIB_LIB_CLUSTER_HPP_

#include <array>
#include <vector>

namespace toppic {

namespace toplib {

namespace lib_cluster {

enum class Linkage { kComplete, kAverage };

// SciPy's scipy.cluster.hierarchy.linkage for the complete and average
// methods: the nearest-neighbour chain algorithm on a condensed distance
// matrix (the row-major upper triangle, i < j), the rows sorted by distance
// (stable) and relabelled as SciPy does, so that equal distances are broken
// the same way. Each row is (left id, right id, distance, size).
std::vector<std::array<double, 4>> linkage(std::vector<double> dist, int n,
                                           Linkage method);

// SciPy's fcluster(Z, t, criterion='distance'): flat clusters whose
// cophenetic distance is at most t, numbered from 1 in SciPy's traversal
// order of the dendrogram.
std::vector<int> fclusterDistance(const std::vector<std::array<double, 4>>& z,
                                  int n, double t);

// The members of each flat cluster (0-based positions, increasing),
// clusters in label order: TopLib's cluster_index_gen on fcluster's labels.
std::vector<std::vector<int>> clusterMembers(const std::vector<int>& labels);

}  // namespace lib_cluster

}  // namespace toplib

}  // namespace toppic

#endif
