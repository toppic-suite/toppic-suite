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


#ifndef ZERO_PTM_FILTER_PROCESSOR_HPP_
#define ZERO_PTM_FILTER_PROCESSOR_HPP_

#include "seq/db_block.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/zero_ptm_filter_mng.hpp"

namespace toppic {

class ZeroPtmFilterProcessor {
 public:
  explicit ZeroPtmFilterProcessor(ZeroPtmFilterMngPtr mng_ptr): mng_ptr_(mng_ptr) {}
  void process();

 private:
  ZeroPtmFilterMngPtr mng_ptr_;

  void processBlock(DbBlockPtr block_ptr);
};

typedef std::shared_ptr<ZeroPtmFilterProcessor> ZeroPtmFilterProcessorPtr;

} /* namespace toppic */

#endif /* ZERO_PTM_FILTER_PROCESSOR_HPP_ */
