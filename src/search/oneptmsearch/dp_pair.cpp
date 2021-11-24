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

#include <limits>

#include "search/oneptmsearch/path_type.hpp"
#include "search/oneptmsearch/dp_pair.hpp"

namespace toppic {

DpPair::DpPair(int x,int y,double pair_score,double diff,
               int order,int n_shift,DiagHeaderPtr header_ptr):Pair(x,y){
  diff_ = diff;
  pair_score_ = pair_score;
  order_= order;
  header_ptr_ = header_ptr;
  for(int i=0;i<n_shift+1;i++){
    prev_pair_ptrs_.push_back(nullptr);
    scores_.push_back(-std::numeric_limits<double>::max());
    types_.push_back(path_type::TYPE_NULL);
  }
}

void DpPair::updateTable(int s,double score,int path_type,DpPairPtr prev_pair_ptr){
  scores_[s] = score;
  types_[s] = path_type;
  prev_pair_ptrs_[s] = prev_pair_ptr;
}

} /* namespace toppic */
