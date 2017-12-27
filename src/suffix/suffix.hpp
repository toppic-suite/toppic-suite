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


#ifndef PROT_SUFFIX_SUFFIX_HPP
#define PROT_SUFFIX_SUFFIX_HPP

#include "node.hpp"

namespace prot {

namespace suffix {

class Suffix {
 public:
  Suffix(NodePtr originNode, int beginIndex, int endIndex):
      origin_node_(originNode),
      bgn_idx_(beginIndex),
      end_idx_(endIndex) {}

  bool isExplicit() {return bgn_idx_ > end_idx_;}

  bool isImplicit() {return end_idx_ >= bgn_idx_;}

  void canonize();

  int getSpan() {return end_idx_ - bgn_idx_;}

  NodePtr getOriginNode() {return origin_node_;}

  int getBeginIndex() {return bgn_idx_;}

  void incBeginIndex() {bgn_idx_++;}

  void changeOriginNode() {origin_node_ = origin_node_->getSuffixNode();}

  int getEndIndex() {return end_idx_;}

  void incEndIndex() {end_idx_++;}

 private:
  NodePtr origin_node_;

  int bgn_idx_;

  int end_idx_;
};

typedef std::shared_ptr<Suffix> SuffixPtr;

}  // namespace suffix

}  // namespace prot
#endif
