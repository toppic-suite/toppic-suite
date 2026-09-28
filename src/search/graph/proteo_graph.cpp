// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "search/graph/proteo_graph.hpp"

#include <algorithm>
#include <cstdlib>
#include <vector>

#include "common/base/mod_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/residue_util.hpp"
#include "common/util/logger.hpp"
#include "seq/alter_type.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/proteoform_factory.hpp"
#include "seq/residue_seq.hpp"

namespace toppic {

ProteoGraph::ProteoGraph(const FastaSubSeqPtr& fasta_seq_ptr,
                         const ModPtrVec& fix_mod_ptr_vec,
                         const MassGraphPtr& graph_ptr, bool is_nme,
                         double convert_ratio, int max_mod_num,
                         int max_ptm_sum_mass, int proteo_graph_gap,
                         int var_ptm_in_gap,
                         const std::vector<char>* dist_cover)
    : is_nme_(is_nme),
      proteo_graph_gap_(proteo_graph_gap),
      var_ptm_in_gap_(var_ptm_in_gap) {
  db_proteo_ptr_ = proteoform_factory::geneDbProteoformPtr(
      fasta_seq_ptr, fix_mod_ptr_vec, fasta_seq_ptr->getSubSeqStart());
  graph_ptr_ = graph_ptr;
  node_num_ = num_vertices(*graph_ptr.get());
  LOG_DEBUG("node num " << node_num_);
  pair_num_ = node_num_ * (proteo_graph_gap_ + 1);
  compSeqMasses(convert_ratio);
  compDistances(max_mod_num, max_ptm_sum_mass, dist_cover);
}

int ProteoGraph::getVecIndex(int v1, int v2) {
  int index = (proteo_graph_gap_ + 1) * v1 + (v2 - v1);
  return index;
}

int ProteoGraph::getSeqMass(int v1, int v2) {
  int index = getVecIndex(v1, v2);
  return seq_masses_[index];
}

void ProteoGraph::compSeqMasses(double convert_ratio) {
  ResSeqPtr res_seq_ptr = db_proteo_ptr_->getResSeqPtr();
  seq_masses_ = std::vector<int>(pair_num_, 0);
  for (int i = 0; i < node_num_; i++) {
    int mass = 0;
    for (int j = i + 1; j < node_num_ && j <= i + proteo_graph_gap_; j++) {
      int cur_mass = std::round(res_seq_ptr->getResiduePtr(j - 1)->getMass() *
                                convert_ratio);
      mass += cur_mass;
      int index = getVecIndex(i, j);
      seq_masses_[index] = mass;
    }
  }
}

// With a cover (SpecGraph::compDistCover), the distance lists keep only the
// distances an alignment with that spectrum graph can pair; the others are
// counted in the lists' all_ fields, from which the alignment takes the
// range of the proteoform's distances.
void ProteoGraph::compDistances(int max_mod_num, int max_ptm_sum_mass,
                                const std::vector<char>* dist_cover) {
  MassGraph* g_p = graph_ptr_.get();
  // The distances of each pair (i, j) with each number of variable
  // modifications, sorted and unique: those of the pair (i, j - 1) plus the
  // mass of an edge from j - 1 to j, collected and then sorted, which gives
  // the set of them.
  std::vector<std::vector<std::vector<int>>> dist_vecs(
      pair_num_, std::vector<std::vector<int>>(max_mod_num + 1));
  // initialize pair (i, i)
  for (int i = 0; i < node_num_; i++) {
    int index = getVecIndex(i, i);
    dist_vecs[index][0].push_back(0);
  }
  for (int i = 0; i < node_num_ - 1; i++) {
    for (int j = i + 1; j < node_num_ && j <= i + proteo_graph_gap_; j++) {
      Vertex v2 = vertex(j, *g_p);
      Vertex pre_v2 = vertex(j - 1, *g_p);
      int index = getVecIndex(i, j);
      int pre_index = getVecIndex(i, j - 1);
      boost::graph_traits<MassGraph>::out_edge_iterator ei, ei_end;
      boost::tie(ei, ei_end) = out_edges(pre_v2, *g_p);
      for (; ei != ei_end; ++ei) {
        if (target(*ei, *g_p) == v2) {
          MassGraph::edge_descriptor e = *ei;
          int d = (*g_p)[e].int_mass_;
          int change = (*g_p)[e].alter_type_;
          bool var_change = change == AlterType::PROTEIN_VARIABLE->getId() ||
                            change == AlterType::VARIABLE->getId();
          for (int k = 0; k < var_ptm_in_gap_ + 1; k++) {
            if (k == max_mod_num && var_change) {
              continue;
            }
            const std::vector<int>& pre_dists = dist_vecs[pre_index][k];
            std::vector<int>& dists = dist_vecs[index][var_change ? k + 1 : k];
            for (size_t t = 0; t < pre_dists.size(); t++) {
              int new_d = d + pre_dists[t];
              if (std::abs(new_d - seq_masses_[index]) <= max_ptm_sum_mass) {
                dists.push_back(new_d);
              }
            }
          }
        }
      }
      for (int k = 0; k < max_mod_num + 1; k++) {
        std::vector<int>& dists = dist_vecs[index][k];
        std::sort(dists.begin(), dists.end());
        dists.erase(std::unique(dists.begin(), dists.end()), dists.end());
      }
    }
  }

  dist_lists_.resize(max_mod_num + 1);
  for (int k = 0; k < max_mod_num + 1; k++) {
    addToDistList(dist_vecs, node_num_, k, dist_lists_[k], proteo_graph_gap_,
                  dist_cover);
  }
}

}  // namespace toppic
