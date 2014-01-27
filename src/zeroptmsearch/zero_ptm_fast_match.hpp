#ifndef ZERO_PTM_FAST_MATCH_HPP_
#define ZERO_PTM_FAST_MATCH_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"

namespace prot {

class ZeroPtmFastMatch {
 public:
  ZeroPtmFastMatch (ProteoformPtr proteoform_ptr, double score, 
                    int begin, int end) {
    proteoform_ptr_ = proteoform_ptr;
    score_ = score;
    begin_ = begin;
    end_ = end;
  }
  double getScore() {return score_;}

  ProteoformPtr getProteoformPtr() {return proteoform_ptr_;}

  int getBegin() {return begin_;}

  int getEnd() {return end_;}

 private:
  ProteoformPtr proteoform_ptr_;
  double score_;
  int begin_;
  int end_;
};

typedef std::shared_ptr<ZeroPtmFastMatch> ZpFastMatchPtr;
typedef std::vector<ZpFastMatchPtr> ZpFastMatchPtrVec;

inline bool compareZeroPtmFastMatchDown(ZpFastMatchPtr m1, 
                                        ZpFastMatchPtr m2) {
  return m1->getScore() > m2->getScore();
}

ZpFastMatchPtrVec zeroPtmFastFilter(SemiAlignTypePtr semi_align_type,
                                    ExtendMsPtr ms_ptr,
                                    ProteoformPtrVec &form_ptr_vec,
                                    int report_num);

}
#endif

