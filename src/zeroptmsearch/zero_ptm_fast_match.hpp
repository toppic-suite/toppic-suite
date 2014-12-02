#ifndef ZERO_PTM_FAST_MATCH_HPP_
#define ZERO_PTM_FAST_MATCH_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"

namespace prot {

class ZeroPtmFastMatch {
 public:
  ZeroPtmFastMatch (ProteoformPtr proteoform_ptr, double score, 
                    int begin, int end);

  double getScore() {return score_;}

  ProteoformPtr getProteoformPtr() {return proteo_ptr_;}

  int getBegin() {return begin_;}

  int getEnd() {return end_;}

 private:
  ProteoformPtr proteo_ptr_;
  double score_;
  int begin_;
  int end_;
};

typedef std::shared_ptr<ZeroPtmFastMatch> ZpFastMatchPtr;
typedef std::vector<ZpFastMatchPtr> ZpFastMatchPtrVec;

inline bool compareZeroPtmFastMatchDown(const ZpFastMatchPtr &a, 
                                        const ZpFastMatchPtr &b) {
  return a->getScore() > b->getScore();
}

ZpFastMatchPtrVec zeroPtmFastFilter(SemiAlignTypePtr semi_align_type_ptr,
                                    const ExtendMsPtrVec &ms_ptr_ptr,
                                    const ProteoformPtrVec &proteo_ptrs,
                                    int report_num);

}
#endif

