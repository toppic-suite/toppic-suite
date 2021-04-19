//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_MS_SPEC_SUPPORT_PEAK_HPP_
#define TOPPIC_MS_SPEC_SUPPORT_PEAK_HPP_

#include "common/base/support_peak_type.hpp"
#include "ms/spec/deconv_peak.hpp"

namespace toppic {

class SupportPeak {
 public:
  SupportPeak(DeconvPeakPtr deconv_peak_ptr, double offset,
              double score, SPTypePtr peak_type_ptr);

  SPTypePtr getPeakTypePtr() {return peak_type_ptr_;}

  double getOffset() {return offset_;}

  double getScore() {return score_;}

  DeconvPeakPtr getDeconvPeakPtr() {return deconv_peak_ptr_;}

 private:
  DeconvPeakPtr deconv_peak_ptr_;

  double offset_;

  double score_;

  SPTypePtr peak_type_ptr_;
};

typedef std::shared_ptr<SupportPeak> SupportPeakPtr;
typedef std::vector<SupportPeakPtr> SupportPeakPtrVec;

} /* namespace toppic */

#endif /* SUPPORT_PEAK_HPP_ */
