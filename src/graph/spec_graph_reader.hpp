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


#ifndef TOPPIC_GRAPH_SPEC_GRAPH_READER_HPP_
#define TOPPIC_GRAPH_SPEC_GRAPH_READER_HPP_

#include "spec/msalign_reader.hpp"
#include "graph/spec_graph.hpp"

namespace toppic {

class SpecGraphReader {
 public:
  SpecGraphReader(const std::string &sp_file_name,
                  int group_sp_num, double convert_ratio,
                  SpParaPtr sp_para_ptr);

  SpecGraphPtrVec getNextSpecGraphPtrVec(int error);

  SpecGraphPtrVec getNextSpecGraphPtrVec(SpectrumSetPtr spec_set_ptr, int error);

 private:
  MsAlignReaderPtr ms_reader_ptr_;
  int group_sp_num_;
  double convert_ratio_;
  SpParaPtr sp_para_ptr_;

  MassGraphPtr getMassGraphPtr(const PrmPeakPtrVec &peak_vec);
};

}  // namespace toppic

#endif /* SPEC_GRAPH_READER_HPP_ */
