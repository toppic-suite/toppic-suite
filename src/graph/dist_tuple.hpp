#ifndef PROT_DIST_TUPLE_HPP_
#define PROT_DIST_TUPLE_HPP_

#include "graph/graph.hpp"

namespace prot {

class DistTuple {
 public:
  DistTuple (MassGraphPtr graph_ptr, int v1, int v2, int mod_num, int dist) {
    graph_ptr_ = graph_ptr;
    v1_ = v1;
    v2_ = v2;
    mod_num_ = mod_num;
    dist_ = dist;
  }

  MassGraphPtr getGraphPtr() {return graph_ptr_;}
  int getFirstVertex() {return v1_;}
  int getSecondVertex() {return v2_;}
  int getModNum() {return mod_num_;}
  int getDist() {return dist_;}

 private:
  MassGraphPtr graph_ptr_;
  int v1_;
  int v2_;
  int mod_num_;
  int dist_;
};

typedef std::shared_ptr<DistTuple> DistTuplePtr;
typedef std::vector<DistTuplePtr> DistTuplePtrVec;
typedef std::vector<DistTuplePtrVec> DistTuplePtrVec2D;

inline bool distUp(const DistTuplePtr &a, const DistTuplePtr &b){
  return a->getDist() < b->getDist();
}

}
#endif
