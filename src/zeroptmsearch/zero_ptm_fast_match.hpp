#ifndef ZERO_PTM_FAST_MATCH_HPP_
#define ZERO_PTM_FAST_MATCH_HPP_

#include "base/proteoform.hpp"

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

}
#endif

