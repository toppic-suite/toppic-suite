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


#ifndef TOPPIC_SEARCH_GRAPH_PROTEO_GRAPH_READER_HPP_
#define TOPPIC_SEARCH_GRAPH_PROTEO_GRAPH_READER_HPP_

#include <string>

#include "search/graph/graph.hpp"
#include "search/graph/proteo_anno.hpp"
#include "search/graph/proteo_graph.hpp"

namespace toppic {

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

 private:
  ModPtrVec fix_mod_ptr_vec_;

  double convert_ratio_;

  int max_mod_num_;

  int max_ptm_sum_mass_;

  ProteoAnnoPtr proteo_anno_ptr_;

  FastaReaderPtr reader_ptr_;

  int proteo_graph_gap_;

  int var_ptm_in_gap_;
};

MassGraphPtr getMassGraphPtr(ProteoAnnoPtr proteo_anno_ptr, double convert_ratio);

} /* namespace toppic */

#endif /* PROTEO_GRAPH_READER_HPP_ */
