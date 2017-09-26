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


#ifndef PROT_GRAPH_ALIGN_HPP_
#define PROT_GRAPH_ALIGN_HPP_

#include "prsm/prsm.hpp"
#include "oneptmsearch/diagonal_header.hpp"
#include "graph/graph.hpp"
#include "graph/proteo_graph.hpp"
#include "graph/spec_graph.hpp"
#include "graphalign/graph_dp_node.hpp"
#include "graphalign/graph_result_node.hpp"
#include "graphalign/graph_align_mng.hpp"

namespace prot {

typedef std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>> ConsistentPairs;

class GraphAlign {
 public:
  GraphAlign(GraphAlignMngPtr mng_ptr, ProteoGraphPtr proteo_graph_ptr,
             SpecGraphPtr spec_graph_ptr);
  void process();

  PrsmPtr geneResult(int s, int m);
  PrsmPtr geneResult(int s);

 private:
  GraphAlignMngPtr mng_ptr_;
  ProteoGraphPtr proteo_graph_ptr_;
  MassGraphPtr pg_;
  int proteo_ver_num_;
  SpecGraphPtr spec_graph_ptr_;
  MassGraphPtr sg_;
  int spec_ver_num_;
  int n_unknown_shift_;
  DistVec spec_dist_;
  DistVec2D dist_vec_;
  ConsistentPairs cons_pairs_;
  GraphDpNodePtrVec2D table_;
  GraphResultNodePtrVec3D result_nodes_;

  GraphResultNodePtrVec2D nodes_2d_;
  DiagonalHeaderPtrVec diag_headers_; 
  DiagonalHeaderPtrVec2D diag_headers_2d_; 

  void getConsistentPairs();
  void addToConsistentPairs(int m, const std::vector<std::pair<int, int>> & sp_pair_ij,
                            const std::vector<std::pair<int, int>> & pg_pair_ij);
  void initTable();

  GraphDpNodePtr compBestVariableNode(int i, int j, int s, int m, int &best_edge_mod_num);
  GraphDpNodePtr compBestShiftNode(int i, int j, int s, int m);
  void updateBestShiftNode(int i, int j, int s, int m);
  void dp();

  GraphResultNodePtrVec backtrace(int s, int m);
  void backtrace();

  void getNodeDiagonals(int s, int m);

  void geneHeaders();

  /*

  //tempary


*/
};

typedef std::shared_ptr<GraphAlign> GraphAlignPtr;

}

#endif

