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
  SupportPeak(const DeconvPeakPtr &deconv_peak_ptr, double offset,
              double score, const SPTypePtr &peak_type_ptr);
  SPTypePtr getPeakTypePtr(){return peak_type_ptr_;}
  double getOffset(){return offset_;}
  double getScore(){return score_;}
  DeconvPeakPtr getDeconvPeakPtr(){return deconv_peak_ptr_;}

 private:
  DeconvPeakPtr deconv_peak_ptr_;
  double offset_;
  double score_;
  SPTypePtr peak_type_ptr_;
};
typedef std::shared_ptr<SupportPeak> SupportPeakPtr;
typedef std::vector<SupportPeakPtr> SupportPeakPtrVec;

} /* namespace prot */

#endif /* SUPPORT_PEAK_HPP_ */
