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

#ifndef TOPPIC_FILTER_ONE_INDEX_ONE_PTM_INDEX_PROCESSOR_HPP_
#define TOPPIC_FILTER_ONE_INDEX_ONE_PTM_INDEX_PROCESSOR_HPP_

#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "seq/db_block.hpp"
namespace toppic{
    
class OnePtmIndexProcessor {
 public:
  explicit OnePtmIndexProcessor(OnePtmFilterMngPtr mng_ptr): mng_ptr_(mng_ptr) {}
  void process();
 
 private:
  OnePtmFilterMngPtr mng_ptr_;
};
typedef std::shared_ptr<OnePtmIndexProcessor> OnePtmIndexProcessorPtr;

}/* namespace toppic */

#endif /* ONE_PTM_INDEX_PROCESSOR_HPP_ */
