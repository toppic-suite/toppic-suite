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

  ProteoformPtr getFormPtr() {return proteoform_ptr_;}

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


inline double compareZeroPtmFastMatchDown(ZpFastMatchPtr m1, ZpFastMatchPtr m2) {
  return m2->getScore() - m1->getScore();
}


ZpFastMatchPtrVec zeroPtmFastFilter(int semi_align_type,
                                    ExtendMsPtr ms_ptr,
                                    ProteoformPtrVec &form_ptr_vec,
                                    int report_num);


ZpFastMatchPtr computeCompMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr);

ZpFastMatchPtr computePrefixMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr);

ZpFastMatchPtr computeSuffixMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr);

ZpFastMatchPtr computeInternalMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr);
/*
 * in the computing of diagonal score in fast filtering, we allow to use n
 * terminal large error tolerance
 */
double compDiagScr(ExtendMsPtr ms_ptr,
                  std::vector<double> &masses, double center);

}
#endif

