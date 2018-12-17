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


#ifndef PROT_SPEC_GRAPH_HPP_
#define PROT_SPEC_GRAPH_HPP_

#include <vector>

#include "spec/prm_peak.hpp"
#include "spec/spectrum_set.hpp"
#include "graph/dist.hpp"
#include "graph/graph.hpp"

namespace toppic {

class SpecGraph {
 public:
  explicit SpecGraph(SpectrumSetPtr spec_set_ptr): spec_set_ptr_(spec_set_ptr) {}

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

}  // namespace toppic

#endif /* SPEC_GRAPH_HPP_ */
