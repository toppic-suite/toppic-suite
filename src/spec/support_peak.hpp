/*
 * support_peak.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#ifndef PROT_SUPPORT_PEAK_HPP_
#define PROT_SUPPORT_PEAK_HPP_

#include "base/support_peak_type.hpp"
#include "spec/deconv_peak.hpp"

namespace prot {

class SupportPeak {
 public:
  SupportPeak(DeconvPeakPtr peak,double offset,double score,
              SPTypePtr peak_type);
  SPTypePtr getPeakType(){return peak_type_;}
  double getOffset(){return offset_;}
  double getScore(){return score_;}
  DeconvPeakPtr getPeak(){return peak_;}

 private:
  DeconvPeakPtr peak_;
  double offset_;
  double score_;
  SPTypePtr peak_type_;
};
typedef std::shared_ptr<SupportPeak> SupportPeakPtr;
typedef std::vector<SupportPeakPtr> SupportPeakPtrVec;

} /* namespace prot */

#endif /* SUPPORT_PEAK_HPP_ */
