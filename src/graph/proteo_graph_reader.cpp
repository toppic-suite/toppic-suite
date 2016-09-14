// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include "base/logger.hpp"
#include "graph/proteo_graph_reader.hpp"

namespace prot {

ProteoGraphReader::ProteoGraphReader(const std::string &db_file_name,
                                     const ModPtrVec &fix_mod_ptr_vec, 
                                     const ProtModPtrVec &prot_mod_ptr_vec,
                                     const ModPtrVec &var_mod_ptr_vec,
                                     double convert_ratio, int max_mod_num,
                                     int max_ptm_sum_mass, int proteo_graph_gap) {
  fix_mod_ptr_vec_ = fix_mod_ptr_vec;
  convert_ratio_ = convert_ratio;
  max_mod_num_ = max_mod_num;
  max_ptm_sum_mass_ = max_ptm_sum_mass;
  reader_ptr_ = FastaReaderPtr(new FastaReader(db_file_name));
  proteo_anno_ptr_ = ProteoAnnoPtr(
      new ProteoAnno(fix_mod_ptr_vec, prot_mod_ptr_vec, var_mod_ptr_vec));
  proteo_graph_gap_ = proteo_graph_gap;
}

MassGraphPtr ProteoGraphReader::getMassGraphPtr() {
  MassGraphPtr graph_ptr = MassGraphPtr(new MassGraph());
  int seq_len  = proteo_anno_ptr_->getLen();
  for (int i = 0; i < seq_len + 1; i++) {
    VertexInfo v(i);
    add_vertex(v, *graph_ptr.get());
  }

  for (int i = 0; i < seq_len; i++) {
    Vertex v1, v2;
    v1 = vertex(i, *graph_ptr.get());
    v2 = vertex(i+1, *graph_ptr.get());
    ResiduePtrVec res_ptr_vec = proteo_anno_ptr_->getResiduePtrVec(i);
    std::vector<int> change_vec = proteo_anno_ptr_->getChangeVec(i);
    for (size_t j=0; j < res_ptr_vec.size(); j++) { 
      EdgeInfo edge_info(res_ptr_vec[j], change_vec[j], convert_ratio_);
      add_edge(v1, v2, edge_info , *graph_ptr.get());
    }
  }
  return graph_ptr;
}

ProteoGraphPtr ProteoGraphReader::getNextProteoGraphPtr() {
  FastaSeqPtr seq_ptr = reader_ptr_->getNextSeq();
  if (seq_ptr.get() == nullptr) {
    return ProteoGraphPtr(nullptr);
  }
  LOG_DEBUG("name " << seq_ptr->getName() << " seq " << seq_ptr->getRawSeq());
  proteo_anno_ptr_->anno(seq_ptr->getRawSeq());
  MassGraphPtr graph_ptr = getMassGraphPtr(); 

  return ProteoGraphPtr(new ProteoGraph(seq_ptr, fix_mod_ptr_vec_, graph_ptr, 
                                        proteo_anno_ptr_->isNme(),
                                        convert_ratio_, max_mod_num_,
                                        max_ptm_sum_mass_, proteo_graph_gap_));
}

ProteoGraphPtr ProteoGraphReader::getProteoGraphPtrBySeq(FastaSeqPtr seq_ptr) {
  proteo_anno_ptr_->anno(seq_ptr->getRawSeq());
  MassGraphPtr graph_ptr = getMassGraphPtr(); 
  return ProteoGraphPtr(new ProteoGraph(seq_ptr, fix_mod_ptr_vec_, graph_ptr, 
                                        proteo_anno_ptr_->isNme(),
                                        convert_ratio_, max_mod_num_,
                                        max_ptm_sum_mass_, proteo_graph_gap_));	
}

}

