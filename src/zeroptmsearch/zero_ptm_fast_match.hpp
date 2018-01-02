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


#ifndef ZERO_PTM_FAST_MATCH_HPP_
#define ZERO_PTM_FAST_MATCH_HPP_

#include <vector>

#include "base/proteoform.hpp"
#include "spec/extend_ms.hpp"

namespace prot {

class ZeroPtmFastMatch;
typedef std::shared_ptr<ZeroPtmFastMatch> ZpFastMatchPtr;
typedef std::vector<ZpFastMatchPtr> ZpFastMatchPtrVec;

class ZeroPtmFastMatch {
 public:
  ZeroPtmFastMatch(ProteoformPtr proteo_ptr, double score, int begin, int end):
      proteo_ptr_(proteo_ptr),
      score_(score),
      begin_(begin),
      end_(end) {}

  double getScore() {return score_;}

  ProteoformPtr getProteoformPtr() {return proteo_ptr_;}

  int getBegin() {return begin_;}

  int getEnd() {return end_;}

  static bool cmpScoreDec(const ZpFastMatchPtr &a, const ZpFastMatchPtr &b) {
    return a->getScore() > b->getScore();
  }

  static ZpFastMatchPtrVec filter(AlignTypePtr align_type_ptr,
                                  const ExtendMsPtrVec &ms_ptr_ptr,
                                  const ProteoformPtrVec &proteo_ptrs,
                                  int report_num, double ppo);

 private:
  ProteoformPtr proteo_ptr_;
  double score_;
  int begin_;
  int end_;
};

}  // namespace prot
#endif

