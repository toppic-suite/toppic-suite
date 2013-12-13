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

inline double compareZeroPtmFastMatchDown(ZeroPtmFastMatch m1, ZeroPtmFastMatch m2) {
  return m2.getScore() - m1.getScore();
}


ZeroPtmFastMatch computeCompMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr);

ZeroPtmFastMatch computePrefixMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr);

/*
 * in the computing of diagonal score in fast filtering, we allow to use n
 * terminal large error tolerance
 */
double compDiagScr(ExtendMsPtr ms_ptr,
                  std::vector<double> &masses, double center);

}
#endif

