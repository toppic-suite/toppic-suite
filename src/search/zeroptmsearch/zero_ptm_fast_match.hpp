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

#ifndef TOPPIC_SEARCH_ZERO_PTM_SEARCH_ZERO_PTM_FAST_MATCH_HPP_
#define TOPPIC_SEARCH_ZERO_PTM_SEARCH_ZERO_PTM_FAST_MATCH_HPP_

#include <vector>

#include "seq/proteoform.hpp"

namespace toppic {

class ZeroPtmFastMatch;
typedef std::shared_ptr<ZeroPtmFastMatch> ZpFastMatchPtr;
typedef std::vector<ZpFastMatchPtr> ZpFastMatchPtrVec;

class ZeroPtmFastMatch {
 public:
  ZeroPtmFastMatch(ProteoformPtr proteo_ptr, double score, int begin, int end);

  double getScore() {return score_;}

  ProteoformPtr getProteoformPtr() {return proteo_ptr_;}

  int getBegin() {return begin_;}

  int getEnd() {return end_;}

  static bool cmpScoreDec(const ZpFastMatchPtr &a, const ZpFastMatchPtr &b) {
    return a->getScore() > b->getScore();
  }

 private:
  ProteoformPtr proteo_ptr_;
  double score_;
  int begin_;
  int end_;
};

}  // namespace toppic
#endif

