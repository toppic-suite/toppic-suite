// Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#ifndef TOPPIC_SEARCH_GRAPH_ALIGN_GRAPH_POST_PROCESSOR_HPP_
#define TOPPIC_SEARCH_GRAPH_ALIGN_GRAPH_POST_PROCESSOR_HPP_

#include <unordered_map>
#include <string>

#include "search/graphalign/graph_align_mng.hpp"

namespace toppic {

class GraphPostProcessor {
 public:
  GraphPostProcessor(GraphAlignMngPtr mng_ptr,
                     const std::string & input_ext,
                     const std::string & output_ext):
      mng_ptr_(mng_ptr),
      input_ext_(input_ext),
      output_ext_(output_ext) {}

  void process();

  void geneMassPtmMap();

 private:
  GraphAlignMngPtr mng_ptr_;

  std::string input_ext_;

  std::string output_ext_;

  std::unordered_map<int, PtmPtrVec> mass_ptm_map_;
};

typedef std::shared_ptr<GraphPostProcessor> GraphPostProcessorPtr;

}  // namespace toppic

#endif
