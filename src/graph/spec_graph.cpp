#include "base/logger.hpp"
#include "graph/graph_util.hpp"
#include "graph/spec_graph.hpp"

namespace prot {

SpecGraph::SpecGraph(SpectrumSetPtr spec_set_ptr) {
  spec_set_ptr_ = spec_set_ptr;
}


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

}

