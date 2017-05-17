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

#include <set>
#include <vector>

#include "base/logger.hpp"
#include "base/change_type.hpp"
#include "base/fasta_reader.hpp"
#include "base/residue_seq.hpp"
#include "base/proteoform_factory.hpp"
#include "graph/graph_util.hpp"
#include "graph/proteo_graph.hpp"

namespace prot {

ProteoGraph::ProteoGraph(FastaSeqPtr fasta_seq_ptr, ModPtrVec fix_mod_ptr_vec,
                         MassGraphPtr graph_ptr, bool is_nme,
                         double convert_ratio, int max_mod_num,
                         int max_ptm_sum_mass, int proteo_graph_gap) {
  db_proteo_ptr_ = ProteoformFactory::geneDbProteoformPtr(fasta_seq_ptr, fix_mod_ptr_vec);
  graph_ptr_ = graph_ptr;
  is_nme_ = is_nme;
  node_num_ = num_vertices(*graph_ptr.get());
  LOG_DEBUG("node num " << node_num_);
  proteo_graph_gap_ = proteo_graph_gap;
  pair_num_ = node_num_ * (proteo_graph_gap_ + 1);
  compSeqMasses(convert_ratio);
  compDistances(max_mod_num, max_ptm_sum_mass);
}

int ProteoGraph::getVecIndex(int v1, int v2) {
  int index =  (proteo_graph_gap_ + 1) * v1 + (v2 - v1);
  return index;
}

int ProteoGraph::getSeqMass(int v1, int v2) {
  int index = getVecIndex(v1, v2);
  return seq_masses_[index];
}

void ProteoGraph::compSeqMasses(double convert_ratio) {
  ResSeqPtr res_seq_ptr = db_proteo_ptr_->getResSeqPtr();
  seq_masses_ = std::vector<int>(pair_num_, 0);
  for (int i = 0; i < node_num_; i ++) {
    int mass = 0;
    for (int j = i + 1; j < node_num_ && j <= i + proteo_graph_gap_; j++) {
      int cur_mass = std::round(res_seq_ptr->getResiduePtr(j-1)->getMass() * convert_ratio);
      mass += cur_mass;
      int index = getVecIndex(i, j);
      seq_masses_[index] = mass;
    }
  }
}

void ProteoGraph::compDistances(int max_mod_num, int max_ptm_sum_mass) {
  MassGraph *g_p = graph_ptr_.get();
  // get mass without ptms

  std::vector<std::vector<std::set<int>>> dist_vecs;
  for (int i = 0; i < pair_num_; i++) {
    std::set<int> empty_set;
    std::vector<std::set<int>> one_pair_vec;
    for (int j = 0; j < max_mod_num + 1; j ++) {
      one_pair_vec.push_back(empty_set);
    }
    dist_vecs.push_back(one_pair_vec);
  }
  // initialize pair (i, i)
  for (int i = 0; i < node_num_; i++) {
    int index = getVecIndex(i, i);
    dist_vecs[index][0].insert(0);
  }
  for (int i = 0; i < node_num_ - 1; i++) {
    for (int j = i + 1; j < node_num_ && j <= i + proteo_graph_gap_; j++) {
      Vertex v2 = vertex(j, *g_p);
      Vertex pre_v2 = vertex(j-1, *g_p);
      int index = getVecIndex(i, j);
      int pre_index = getVecIndex(i, j-1);
      boost::graph_traits<MassGraph>::out_edge_iterator ei, ei_end;
      boost::tie(ei, ei_end) = out_edges(pre_v2, *g_p);
      for ( ; ei != ei_end; ++ei) {
        if (target(*ei, *g_p) == v2) {
          MassGraph::edge_descriptor e = *ei;
          int d =(*g_p)[e].int_mass_;
          int change = (*g_p)[e].change_type_;
          for (int k = 0; k < max_mod_num + 1; k++) {
            if (k == max_mod_num &&
                (change == ChangeType::PROTEIN_VARIABLE->getId()
                 || change == ChangeType::VARIABLE->getId())) {
              continue;
            }
            for (std::set<int>::iterator it=dist_vecs[pre_index][k].begin();
                 it != dist_vecs[pre_index][k].end(); it++) {
              int new_d = d + *it;
              if (std::abs(new_d - seq_masses_[index]) <= max_ptm_sum_mass) {
                if (change == ChangeType::PROTEIN_VARIABLE->getId()
                    || change == ChangeType::VARIABLE->getId()) {
                  dist_vecs[index][k+1].insert(new_d);
                } else {
                  dist_vecs[index][k].insert(new_d);
                }
              }
            }
          }
        }
      }
    }
  }

  dist_vec_.reserve(max_mod_num);
  DistVec tmp;
  for (int k = 0; k < max_mod_num + 1; k++) {
    dist_vec_.push_back(tmp);
  }

  for (int k = 0; k < max_mod_num + 1; k++) {
    addToDistVec(graph_ptr_, dist_vecs, node_num_, k, dist_vec_[k], proteo_graph_gap_);
  }
}

}  // namespace prot

