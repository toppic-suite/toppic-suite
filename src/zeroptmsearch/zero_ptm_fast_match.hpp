#ifndef ZERO_PTM_FAST_MATCH_HPP_
#define ZERO_PTM_FAST_MATCH_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"

namespace prot {

class ZeroPtmFastMatch {
 public:
  ZeroPtmFastMatch (ProteoformPtr proteoform_ptr, double score) {
    proteoform_ptr_ = proteoform_ptr;
    score_ = score_;
  }
  double getScore() {return score_;}

  ProteoformPtr getFormPtr() {return proteoform_ptr_;}

 private:
  ProteoformPtr proteoform_ptr_;
  double score_;
};

ZeroPtmFastMatch computeCompMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr);

/*
 * in the computing of diagonal score in fast filtering, we allow to use n
 * terminal large error tolerance
 */
double compDiagScr(ExtendMsPtr ms_ptr,
                  std::vector<double> &masses, double center);

}
#endif

