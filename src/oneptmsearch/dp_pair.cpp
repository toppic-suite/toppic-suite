//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "oneptmsearch/dp_pair.hpp"

namespace prot {
DPPair::DPPair(int x,int y,double pair_score,double diff,
               int order,int n_shift,DiagonalHeaderPtr header_ptr):Pair(x,y){
  diff_ = diff;
  pair_score_ = pair_score;
  order_= order;
  header_ptr_ = header_ptr;
  for(int i=0;i<n_shift+1;i++){
    prev_pair_ptrs_.push_back(nullptr);
    scores_.push_back(-std::numeric_limits<double>::max());
    types_.push_back(PATH_TYPE_NULL);
  }
}

void DPPair::updateTable(int s,double score,int path_type,DPPairPtr prev_pair_ptr){
    scores_[s] = score;
    types_[s] = path_type;
    prev_pair_ptrs_[s] = prev_pair_ptr;
}

} /* namespace prot */
