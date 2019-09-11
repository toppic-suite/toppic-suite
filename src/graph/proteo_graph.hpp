//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#ifndef TOPPIC_GRAPH_PROTEO_GRAPH_HPP_
#define TOPPIC_GRAPH_PROTEO_GRAPH_HPP_

#include "seq/fasta_seq.hpp"
#include "seq/proteoform.hpp"
#include "graph/dist.hpp"
#include "graph/graph.hpp"

namespace toppic {

class ProteoGraph {
 public:
  ProteoGraph(FastaSubSeqPtr seq_ptr, ModPtrVec fix_mod_ptr_vec,
              MassGraphPtr graph_ptr, bool is_nme, 
              double convert_ratio, int max_mod_num,
              int max_ptm_sum_mass, int proteo_graph_gap,
              int var_ptm_in_gap);

  int getVecIndex(int v1, int v2);

  int getSeqMass(int v1, int v2);

  ProteoformPtr getProteoformPtr() {return db_proteo_ptr_;}

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

} /* namespace toppic */

#endif /* PROTEO_GRAPH_HPP_ */
