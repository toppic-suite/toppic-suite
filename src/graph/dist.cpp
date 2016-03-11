
#include "dist.hpp"

namespace prot {

int getVecIndex(int v1, int v2, int node_num) {
  int index =  (node_num + node_num - v1 + 1) * v1 /2  + (v2 - v1);
  return index;
}

void addToDistVec(MassGraphPtr graph_ptr, const std::vector<std::vector<std::set<int>>> & dist_vecs,
                  int node_num, int mod_num, DistVec & dist_vec, int gap) {
  std::set<Dist> dist_set;
  for (int i = 0; i < node_num - 1; i++) {
    for (int j = i + 1; j < node_num && j <= i + gap; j++) {
      int index = getVecIndex(i, j, node_num);
      for (std::set<int>::iterator it=dist_vecs[index][mod_num].begin(); 
           it!=dist_vecs[index][mod_num].end(); it++) {
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

}
