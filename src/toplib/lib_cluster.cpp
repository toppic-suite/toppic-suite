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

#include "toplib/lib_cluster.hpp"

#include <algorithm>
#include <limits>
#include <numeric>

namespace toppic {

namespace toplib {

namespace lib_cluster {

namespace {

// SciPy's condensed_index
inline size_t condensedIndex(int n, int i, int j) {
  if (i < j) {
    return static_cast<size_t>(n) * i - (static_cast<size_t>(i) * (i + 1) / 2) +
           (j - i - 1);
  }
  return static_cast<size_t>(n) * j - (static_cast<size_t>(j) * (j + 1) / 2) +
         (i - j - 1);
}

// SciPy's _hierarchy_distance_update
inline double newDist(double d_xi, double d_yi, int size_x, int size_y,
                      Linkage method) {
  if (method == Linkage::kComplete) {
    return std::max(d_xi, d_yi);
  }
  return (size_x * d_xi + size_y * d_yi) / (size_x + size_y);
}

// SciPy's LinkageUnionFind
class UnionFind {
 public:
  explicit UnionFind(int n)
      : parent_(2 * n - 1), size_(2 * n - 1, 1), next_label_(n) {
    std::iota(parent_.begin(), parent_.end(), 0);
  }

  int merge(int x, int y) {
    parent_[x] = next_label_;
    parent_[y] = next_label_;
    int size = size_[x] + size_[y];
    size_[next_label_] = size;
    next_label_++;
    return size;
  }

  int find(int x) {
    int p = x;
    while (parent_[x] != x) {
      x = parent_[x];
    }
    while (parent_[p] != x) {
      int next = parent_[p];
      parent_[p] = x;
      p = next;
    }
    return x;
  }

 private:
  std::vector<int> parent_;
  std::vector<int> size_;
  int next_label_;
};

}  // namespace

std::vector<std::array<double, 4>> linkage(std::vector<double> dist, int n,
                                           Linkage method) {
  std::vector<std::array<double, 4>> z(n - 1);
  std::vector<int> size(n, 1);
  std::vector<int> chain(n);
  int chain_length = 0;
  for (int k = 0; k < n - 1; k++) {
    if (chain_length == 0) {
      chain_length = 1;
      for (int i = 0; i < n; i++) {
        if (size[i] > 0) {
          chain[0] = i;
          break;
        }
      }
    }
    int x = 0;
    int y = 0;
    double current_min = 0;
    // go through the chain of neighbours until two mutual neighbours
    while (true) {
      x = chain[chain_length - 1];
      // prefer the previous element of the chain as the minimum
      if (chain_length > 1) {
        y = chain[chain_length - 2];
        current_min = dist[condensedIndex(n, x, y)];
      } else {
        current_min = std::numeric_limits<double>::infinity();
      }
      for (int i = 0; i < n; i++) {
        if (size[i] == 0 || x == i) {
          continue;
        }
        double d = dist[condensedIndex(n, x, i)];
        if (d < current_min) {
          current_min = d;
          y = i;
        }
      }
      if (chain_length > 1 && y == chain[chain_length - 2]) {
        break;
      }
      chain[chain_length] = y;
      chain_length++;
    }
    chain_length -= 2;
    if (x > y) {
      std::swap(x, y);
    }
    int nx = size[x];
    int ny = size[y];
    z[k] = {static_cast<double>(x), static_cast<double>(y), current_min,
            static_cast<double>(nx + ny)};
    size[x] = 0;
    size[y] = nx + ny;
    for (int i = 0; i < n; i++) {
      int ni = size[i];
      if (ni == 0 || i == y) {
        continue;
      }
      dist[condensedIndex(n, i, y)] =
          newDist(dist[condensedIndex(n, i, x)], dist[condensedIndex(n, i, y)],
                  nx, ny, method);
    }
  }
  // sort by distance (numpy mergesort, stable) and relabel (SciPy's label)
  std::stable_sort(z.begin(), z.end(),
                   [](const std::array<double, 4>& a,
                      const std::array<double, 4>& b) { return a[2] < b[2]; });
  UnionFind uf(n);
  for (int i = 0; i < n - 1; i++) {
    int xr = uf.find(static_cast<int>(z[i][0]));
    int yr = uf.find(static_cast<int>(z[i][1]));
    if (xr < yr) {
      z[i][0] = xr;
      z[i][1] = yr;
    } else {
      z[i][0] = yr;
      z[i][1] = xr;
    }
    z[i][3] = uf.merge(xr, yr);
  }
  return z;
}

std::vector<int> fclusterDistance(const std::vector<std::array<double, 4>>& z,
                                  int n, double t) {
  std::vector<int> labels(n, 0);
  if (n == 1) {
    labels[0] = 1;
    return labels;
  }
  // SciPy's get_max_dist_for_each_cluster: the largest merge distance in
  // each subtree
  std::vector<double> max_dist(n - 1, 0);
  {
    std::vector<int> stack(n);
    std::vector<bool> visited(2 * n - 1, false);
    int k = 0;
    stack[0] = 2 * n - 2;
    while (k >= 0) {
      int root = stack[k] - n;
      int lc = static_cast<int>(z[root][0]);
      int rc = static_cast<int>(z[root][1]);
      if (lc >= n && !visited[lc]) {
        visited[lc] = true;
        stack[++k] = lc;
        continue;
      }
      if (rc >= n && !visited[rc]) {
        visited[rc] = true;
        stack[++k] = rc;
        continue;
      }
      double m = z[root][2];
      if (lc >= n) {
        m = std::max(m, max_dist[lc - n]);
      }
      if (rc >= n) {
        m = std::max(m, max_dist[rc - n]);
      }
      max_dist[root] = m;
      k--;
    }
  }
  // SciPy's cluster_monocrit
  std::vector<int> curr_node(n);
  std::vector<bool> visited(2 * n - 1, false);
  int k = 0;
  curr_node[0] = 2 * n - 2;
  int n_cluster = 0;
  int cluster_leader = -1;
  while (k >= 0) {
    int root = curr_node[k] - n;
    int lc = static_cast<int>(z[root][0]);
    int rc = static_cast<int>(z[root][1]);
    if (cluster_leader == -1 && max_dist[root] <= t) {
      cluster_leader = root;
      n_cluster++;
    }
    if (lc >= n && !visited[lc]) {
      visited[lc] = true;
      curr_node[++k] = lc;
      continue;
    }
    if (rc >= n && !visited[rc]) {
      visited[rc] = true;
      curr_node[++k] = rc;
      continue;
    }
    if (lc < n) {
      if (cluster_leader == -1) {
        n_cluster++;
      }
      labels[lc] = n_cluster;
    }
    if (rc < n) {
      if (cluster_leader == -1) {
        n_cluster++;
      }
      labels[rc] = n_cluster;
    }
    if (cluster_leader == root) {
      cluster_leader = -1;
    }
    k--;
  }
  return labels;
}

std::vector<std::vector<int>> clusterMembers(const std::vector<int>& labels) {
  int max_label = 0;
  for (int l : labels) {
    max_label = std::max(max_label, l);
  }
  std::vector<std::vector<int>> members(max_label);
  for (size_t i = 0; i < labels.size(); i++) {
    members[labels[i] - 1].push_back(static_cast<int>(i));
  }
  // cluster_index_gen drops empty labels
  members.erase(std::remove_if(members.begin(), members.end(),
                               [](const std::vector<int>& m) { return m.empty(); }),
                members.end());
  return members;
}

}  // namespace lib_cluster

}  // namespace toplib

}  // namespace toppic
