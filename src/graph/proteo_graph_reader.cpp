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


#include <string>
#include <vector>

#include "common/util/logger.hpp"
#include "graph/proteo_graph_reader.hpp"

namespace toppic {

ProteoGraphReader::ProteoGraphReader(const std::string &db_file_name,
                                     const ModPtrVec &fix_mod_ptr_vec,
                                     const ProtModPtrVec &prot_mod_ptr_vec,
                                     const ModPtrVec &var_mod_ptr_vec,
                                     double convert_ratio,
                                     int max_mod_num,
                                     int max_ptm_sum_mass,
                                     int proteo_graph_gap,
                                     int var_ptm_in_gap):
    fix_mod_ptr_vec_(fix_mod_ptr_vec),
    convert_ratio_(convert_ratio),
    max_mod_num_(max_mod_num),
    max_ptm_sum_mass_(max_ptm_sum_mass),
    proteo_graph_gap_(proteo_graph_gap),
    var_ptm_in_gap_(var_ptm_in_gap) {
      reader_ptr_ = std::make_shared<FastaReader>(db_file_name);
      proteo_anno_ptr_
          = std::make_shared<ProteoAnno>(fix_mod_ptr_vec, prot_mod_ptr_vec, var_mod_ptr_vec);
    }

MassGraphPtr getMassGraphPtr(ProteoAnnoPtr proteo_anno_ptr, double convert_ratio) {
  MassGraphPtr graph_ptr = std::make_shared<MassGraph>();
  int seq_len = proteo_anno_ptr->getLen();
  for (int i = 0; i < seq_len + 1; i++) {
    VertexInfo v(i);
    add_vertex(v, *graph_ptr.get());
  }

  for (int i = 0; i < seq_len; i++) {
    Vertex v1, v2;
    v1 = vertex(i, *graph_ptr.get());
    v2 = vertex(i + 1, *graph_ptr.get());
    ResiduePtrVec res_ptr_vec = proteo_anno_ptr->getResiduePtrVec(i);
    std::vector<int> change_vec = proteo_anno_ptr->getChangeVec(i);
    if (std::find(change_vec.begin(), change_vec.end(), MassShiftType::FIXED->getId()) != change_vec.end()) {
      for (size_t j = 0; j < res_ptr_vec.size(); j++) {
        if (change_vec[j] == MassShiftType::FIXED->getId()) {
          EdgeInfo edge_info(res_ptr_vec[j], change_vec[j], convert_ratio);
          add_edge(v1, v2, edge_info , *graph_ptr.get());
        }
      }
    } else {
      for (size_t j = 0; j < res_ptr_vec.size(); j++) {
        EdgeInfo edge_info(res_ptr_vec[j], change_vec[j], convert_ratio);
        add_edge(v1, v2, edge_info , *graph_ptr.get());
      }
    }
  }
  return graph_ptr;
}

}  // namespace toppic

