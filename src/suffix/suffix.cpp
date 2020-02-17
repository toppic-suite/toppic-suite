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


#include "edge.hpp"
#include "suffix.hpp"

namespace toppic {

namespace suffix {

// match the suffix along the tree from Suffix.originNode
void Suffix::canonize() {
  if (!isExplicit()) {
    EdgePtr edge = origin_node_->findEdge(origin_node_->charAt(bgn_idx_));

    int edgeSpan = edge->getSpan();
    while (edgeSpan <= getSpan()) {
      bgn_idx_ += edgeSpan + 1;
      origin_node_ = edge->getEndNode();
      if (bgn_idx_ <= end_idx_) {
        edge = edge->getEndNode()->findEdge(origin_node_->charAt(bgn_idx_));
        edgeSpan = edge->getSpan();
      }
    }
  }
}

}  // namespace suffix

}  // namespace toppic

