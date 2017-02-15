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
