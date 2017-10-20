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


#ifndef PROT_PROTEO_GRAPH_READER_HPP_
#define PROT_PROTEO_GRAPH_READER_HPP_

#include "graph/graph.hpp"
#include "graph/proteo_anno.hpp"
#include "graph/proteo_graph.hpp"

namespace prot {

class ProteoGraphReader {
 public:
  ProteoGraphReader(const std::string &db_file_name,
                    const ModPtrVec &fix_mod_ptr_vec, 
                    const ProtModPtrVec &prot_mod_ptr_vec,
                    const ModPtrVec &var_mod_ptr_vec,
                    double convert_ratio,
                    int max_mod_num,
                    int max_ptm_sum_mass,
                    int proteo_graph_gap,
                    int var_ptm_in_gap);

  ProteoGraphPtr getNextProteoGraphPtr();

  ProteoGraphPtr getProteoGraphPtrBySeq(FastaSeqPtr seq_ptr);

 private:
  ModPtrVec  fix_mod_ptr_vec_; 
  double convert_ratio_;
  int max_mod_num_;
  int max_ptm_sum_mass_;
  ProteoAnnoPtr proteo_anno_ptr_;
  FastaReaderPtr reader_ptr_;
  int proteo_graph_gap_;
  int var_ptm_in_gap_;

  MassGraphPtr getMassGraphPtr();
};

} /* namespace prot */

#endif /* PROTEO_GRAPH_READER_HPP_ */
