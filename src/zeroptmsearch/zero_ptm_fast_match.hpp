#ifndef ZERO_PTM_FAST_MATCH_HPP_
#define ZERO_PTM_FAST_MATCH_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_ms.hpp"

namespace prot {

class ZeroPtmFastMatch;
typedef std::shared_ptr<ZeroPtmFastMatch> ZpFastMatchPtr;
typedef std::vector<ZpFastMatchPtr> ZpFastMatchPtrVec;

class ZeroPtmFastMatch {
 public:
  ZeroPtmFastMatch (ProteoformPtr proteoform_ptr, double score, 
                    int begin, int end);

  double getScore() {return score_;}

  ProteoformPtr getProteoformPtr() {return proteo_ptr_;}

  int getBegin() {return begin_;}

  int getEnd() {return end_;}

  static bool cmpScoreDec(const ZpFastMatchPtr &a, 
                          const ZpFastMatchPtr &b) {
    return a->getScore() > b->getScore();
  }

  static ZpFastMatchPtrVec filter(AlignTypePtr align_type_ptr,
                                  const ExtendMsPtrVec &ms_ptr_ptr,
                                  const ProteoformPtrVec &proteo_ptrs,
                                  int report_num);
 private:
  ProteoformPtr proteo_ptr_;
  double score_;
  int begin_;
  int end_;
};

}
#endif

