#ifndef PROT_DIST_HPP_
#define PROT_DIST_HPP_

#include "graph/graph.hpp"

namespace prot {

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

inline bool distVecUp(const Dist &a, const Dist &b){
  return a.dist_ < b.dist_;
}

void addToDistVec(MassGraphPtr graph_ptr,
                  const std::vector<std::vector<std::set<int>>> & dist_vecs,
                  int node_num, int mod_num, DistVec & dist_vec,
                  int gap);

}

#endif
