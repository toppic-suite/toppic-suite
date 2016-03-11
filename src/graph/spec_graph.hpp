#ifndef PROT_SPEC_GRAPH_HPP_
#define PROT_SPEC_GRAPH_HPP_

#include "spec/prm_peak.hpp"
#include "spec/spectrum_set.hpp"
#include "graph/dist.hpp"
#include "graph/graph.hpp"

namespace prot {

class SpecGraph {
 public:
  SpecGraph(SpectrumSetPtr spec_set_ptr);
  SpecGraph(SpectrumSetPtr spec_set_ptr, PrmPeakPtrVec peak_vec, 
            MassGraphPtr mass_graph_ptr, double convert_ratio);
  SpectrumSetPtr getSpectrumSetPtr() {return spec_set_ptr_;}
  MassGraphPtr getMassGraphPtr() {return graph_ptr_;}

  DistVec getDistVec() {return dist_;}
  PrmPeakPtr getPrmPeakPtr(int i) {return peak_vec_[i];}
  const PrmPeakPtrVec& getPrmPeakPtrVec() {return peak_vec_;}

  int getPeakDist(int v1, int v2);

 private:
  SpectrumSetPtr spec_set_ptr_;
  int node_num_;
  int pair_num_;
  std::vector<int> peak_dists_;
  MassGraphPtr graph_ptr_;
  PrmPeakPtrVec peak_vec_;

  DistVec dist_;

  int getVecIndex(int v1, int v2);


  void compSpecDistances(double convert_ratio);
};

typedef std::shared_ptr<SpecGraph> SpecGraphPtr;
typedef std::vector<SpecGraphPtr> SpecGraphPtrVec;

} /* namespace prot */

#endif /* SPEC_GRAPH_HPP_ */
