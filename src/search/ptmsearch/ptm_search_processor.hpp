//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_PTM_SEARCH_PTM_SEARCH_PROCESSOR_HPP_
#define TOPPIC_PTM_SEARCH_PTM_SEARCH_PROCESSOR_HPP_

#include "seq/proteoform.hpp"
#include "prsm/simple_prsm.hpp"
#include "search/oneptmsearch/ptm_search_mng.hpp"

namespace toppic {

class PtmSearchProcessor {
 public:
  PtmSearchProcessor(PtmSearchMngPtr mng_ptr);

  void process();

 private:
  PtmSearchMngPtr mng_ptr_;

  ProteoformPtrVec proteo_ptrs_;
  ProteoformPtrVec2D mod_proteo_2d_ptrs_;
  SimplePrsmPtrVec simple_prsm_ptrs_;

};

typedef std::shared_ptr<PtmSearchProcessor> PtmSearchProcessorPtr;

} /* namespace toppic */

#endif /* PTM_SEARCH_PROCESSOR_HPP_ */
