/*
 * support_peak.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include <spec/support_peak.hpp>

namespace prot {
SupportPeak::SupportPeak(const DeconvPeakPtr &peak,double offset,
                         double score, const SPTypePtr &peak_type){
  peak_=peak;
  offset_=offset;
  score_=score;
  peak_type_ = peak_type;
}

} /* namespace prot */
