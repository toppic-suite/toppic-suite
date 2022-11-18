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

#ifndef TOPPIC_SEARCH_VAR_PTM_SEARCH_VAR_PTM_SEARCH_PROCESSOR_HPP_
#define TOPPIC_SEARCH_VAR_PTM_SEARCH_VAR_PTM_SEARCH_PROCESSOR_HPP_

#include "ms/spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "search/varptmsearch/var_ptm_search_mng.hpp"

namespace toppic {

class VarPtmSearchProcessor {
 public:
  explicit VarPtmSearchProcessor(VarPtmSearchMngPtr mng_ptr): mng_ptr_(mng_ptr) {}
  void process();

 private:
  PrsmPtrVec varPtmSearchOneSpec(SpectrumSetPtr spec_set_ptr,
                                 const SimplePrsmPtrVec &simple_prsm_ptr_vec,
                                 FastaIndexReaderPtr reader_ptr,
                                 VarPtmSearchMngPtr mng_ptr,
                                 ProteoformTypePtr type_ptr);

  VarPtmSearchMngPtr mng_ptr_;
};

typedef std::shared_ptr<VarPtmSearchProcessor> VarPtmSearchProcessorPtr;

}  // namespace toppic

#endif
