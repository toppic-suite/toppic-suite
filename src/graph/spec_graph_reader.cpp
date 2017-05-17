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


#include "graph/spec_graph_reader.hpp"

namespace prot {

SpecGraphReader::SpecGraphReader(const std::string &sp_file_name,
                                 int group_sp_num,
                                 double convert_ratio,
                                 SpParaPtr sp_para_ptr) {
  ms_reader_ptr_ = std::make_shared<MsAlignReader>(sp_file_name, group_sp_num,
                                                   sp_para_ptr->getActivationPtr());
  group_sp_num_ = group_sp_num;
  convert_ratio_ = convert_ratio;
  sp_para_ptr_ = sp_para_ptr;
}

MassGraphPtr SpecGraphReader::getMassGraphPtr(const PrmPeakPtrVec &peak_vec) {

  LOG_DEBUG("start mass graph");
  MassGraphPtr graph_ptr = std::make_shared<MassGraph>();

  // add mass 0/start nod
  VertexInfo v(0);
  add_vertex(v, *graph_ptr.get());

  for (size_t i = 1; i < peak_vec.size(); i++) {
    // add a new node for the prm
    VertexInfo cur_v(i);
    add_vertex(cur_v, *graph_ptr.get());

    Vertex v1, v2;
    v1 = vertex(i-1, *graph_ptr.get());
    v2 = vertex(i, *graph_ptr.get());

    double dist = peak_vec[i]->getMonoMass() - peak_vec[i-1]->getMonoMass();

    EdgeInfo edge_info(dist, convert_ratio_);
    add_edge(v1, v2, edge_info , *graph_ptr.get());
  }

  return graph_ptr;
}

SpecGraphPtrVec SpecGraphReader::getNextSpecGraphPtrVec(int error) {
  SpectrumSetPtr spec_set_ptr = ms_reader_ptr_->getNextSpectrumSet(sp_para_ptr_);
  return getNextSpecGraphPtrVec(spec_set_ptr, error);
}

SpecGraphPtrVec SpecGraphReader::getNextSpecGraphPtrVec(SpectrumSetPtr spec_set_ptr, int error) {
  LOG_DEBUG("get spec set ");
  SpecGraphPtrVec graph_ptr_vec;
  if (spec_set_ptr  == nullptr) {
    return graph_ptr_vec;
  }
  std::vector<double> prec_errors;
  prec_errors.push_back(0);
  for (int i = 1; i <= error; i++) {
    prec_errors.push_back(- i * MassConstant::getIsotopeMass());
    prec_errors.push_back(i * MassConstant::getIsotopeMass());
  }

  DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
  //LOG_DEBUG("deconv ms size " << deconv_ms_ptr_vec.size());
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  LOG_DEBUG("prec_mono_mass  " << prec_mono_mass);
  if (spec_set_ptr->isValid()) {
    LOG_DEBUG("valid");
    for (size_t i = 0; i < prec_errors.size(); i++) {
      SpectrumSetPtr adjusted_spec_set_ptr
          = std::make_shared<SpectrumSet>(deconv_ms_ptr_vec, sp_para_ptr_, prec_mono_mass + prec_errors[i]);
      PrmMsPtrVec ms_six_vec = adjusted_spec_set_ptr->getMsSixPtrVec();
      PrmPeakPtrVec peak_vec = PrmMs::getPrmPeakPtrs(ms_six_vec, sp_para_ptr_->getPeakTolerancePtr());
      MassGraphPtr graph_ptr = getMassGraphPtr(peak_vec); 
      LOG_DEBUG("graph complete");
      SpecGraphPtr spec_graph_ptr
          = std::make_shared<SpecGraph>(adjusted_spec_set_ptr, peak_vec, graph_ptr, convert_ratio_);
      graph_ptr_vec.push_back(spec_graph_ptr);
    }
  } else {
    LOG_DEBUG("no valid");
    SpecGraphPtr spec_graph_ptr = std::make_shared<SpecGraph>(spec_set_ptr);
    graph_ptr_vec.push_back(spec_graph_ptr);
  }
  LOG_DEBUG("set geneneted");
  return graph_ptr_vec;
}

}

