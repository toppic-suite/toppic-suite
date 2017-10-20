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


#ifndef PROT_PROTEO_GRAPH_HPP_
#define PROT_PROTEO_GRAPH_HPP_

#include "base/fasta_seq.hpp"
#include "base/proteoform.hpp"
#include "graph/dist.hpp"
#include "graph/graph.hpp"

namespace prot {

class ProteoGraph {
 public:
  ProteoGraph(FastaSeqPtr seq_ptr, ModPtrVec fix_mod_ptr_vec,
              MassGraphPtr graph_ptr, bool is_nme, 
              double convert_ratio, int max_mod_num,
              int max_ptm_sum_mass, int proteo_graph_gap,
              int var_ptm_in_gap);

  int getVecIndex(int v1, int v2);
  int getSeqMass(int v1, int v2);

  ProteoformPtr getProteoformPtr() {return db_proteo_ptr_;}
  FastaSeqPtr getFastaSeqPtr() {return db_proteo_ptr_->getFastaSeqPtr();}
  MassGraphPtr getMassGraphPtr() {return graph_ptr_;}
  bool isNme() {return is_nme_;}

  const DistVec2D& getDistVec2D() {return dist_vec_;}

 private:
  ProteoformPtr db_proteo_ptr_;
  int node_num_;
  int pair_num_;
  std::vector<int> seq_masses_;
  MassGraphPtr graph_ptr_;
  bool is_nme_;

  int proteo_graph_gap_;
  int var_ptm_in_gap_;

  DistVec2D dist_vec_;

  void compSeqMasses(double convert_ratio);
  void compDistances(int max_mod_num, int max_ptm_sum_mass);
};

typedef std::shared_ptr<ProteoGraph> ProteoGraphPtr;
typedef std::vector<ProteoGraphPtr> ProteoGraphPtrVec;

} /* namespace prot */

#endif /* PROTEO_GRAPH_HPP_ */
