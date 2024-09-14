//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include <algorithm>

#include "ms/spec/peak.hpp"
#include "ms/spec/ms.hpp"
#include "ms/spec/deconv_ms_util.hpp"

namespace toppic {

namespace deconv_ms_util {

DeconvMsPtrVec getRefineMsPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                 double new_prec_mass) {
  DeconvMsPtrVec result_ptrs;
  for (std::size_t m = 0; m < deconv_ms_ptr_vec.size(); m++) {
    DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[m];
    MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
    MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mass);
    std::vector<DeconvPeakPtr> peak_ptr_list;
    for (std::size_t p = 0; p < deconv_ms_ptr->size(); p++) {
      DeconvPeakPtr ori_peak_ptr = deconv_ms_ptr->getPeakPtr(p);
      // * is a dereference operator
      DeconvPeakPtr new_peak_ptr = std::make_shared<DeconvPeak>(*ori_peak_ptr.get());
      peak_ptr_list.push_back(new_peak_ptr);
    }
    DeconvMsPtr ms_ptr = std::make_shared<Ms<DeconvPeakPtr> >(header_ptr, peak_ptr_list);
    result_ptrs.push_back(ms_ptr);
  }
  return result_ptrs;
}

void keepTopPeaks(DeconvMsPtrVec &deconv_ms_ptr_vec, size_t peak_num) {
  for (std::size_t m = 0; m < deconv_ms_ptr_vec.size(); m++) {
    DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[m];
    std::vector<DeconvPeakPtr> peak_ptr_list;
    for (std::size_t p = 0; p < deconv_ms_ptr->size(); p++) {
      if (p >= peak_num) {
        break;
      }
      DeconvPeakPtr peak_ptr = deconv_ms_ptr->getPeakPtr(p);
      peak_ptr_list.push_back(peak_ptr);
    }
    deconv_ms_ptr->setPeakPtrVec(peak_ptr_list);
  }
}

void log2Transform(DeconvMsPtrVec &deconv_ms_ptr_vec) { 
  for (std::size_t m = 0; m < deconv_ms_ptr_vec.size(); m++) {
    DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[m];
    for (std::size_t p = 0; p < deconv_ms_ptr->size(); p++) {
      DeconvPeakPtr peak_ptr = deconv_ms_ptr->getPeakPtr(p);
      peak_ptr->log2Transform();
    }
  }
}

void vectorNorm(DeconvMsPtrVec &deconv_ms_ptr_vec) { 
  for (std::size_t m = 0; m < deconv_ms_ptr_vec.size(); m++) {
    DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[m];
    double sum_square = 0;
    for (std::size_t p = 0; p < deconv_ms_ptr->size(); p++) {
      double inte = deconv_ms_ptr->getPeakPtr(p)->getIntensity();
      sum_square = sum_square + (inte * inte);
    }
    double vector_len = std::sqrt(sum_square);
    if (vector_len > 0) {
      for (std::size_t p = 0; p < deconv_ms_ptr->size(); p++) {
        double inte = deconv_ms_ptr->getPeakPtr(p)->getIntensity();
        inte = inte /vector_len; 
        deconv_ms_ptr->getPeakPtr(p)->setIntensity(inte);
      }
    }
  }
}

double compDotProd(DeconvMsPtr ms1, DeconvMsPtr ms2, double ppo) {
  DeconvPeakPtrVec peak_list_1 = ms1->getPeakPtrVec();
  DeconvPeakPtrVec peak_list_2 = ms2->getPeakPtrVec();
  std::sort(peak_list_1.begin(), peak_list_1.end(), Peak::cmpPosInc);
  std::sort(peak_list_2.begin(), peak_list_2.end(), Peak::cmpPosInc);
  int i = 0; 
  int j = 0;
  int len_1 = peak_list_1.size();
  int len_2 = peak_list_2.size();
  double dot_prod = 0;
  while (i < len_1 && j < len_2) { 
    DeconvPeakPtr peak_1 = peak_list_1[i];
    DeconvPeakPtr peak_2 = peak_list_2[j];
    double tol = peak_1->getPosition() * ppo;
    if (peak_1->getCharge() == peak_2->getCharge() && 
        std::abs(peak_1->getPosition() - peak_2->getPosition()) <= tol) {
      dot_prod = dot_prod + peak_1->getIntensity() * peak_2->getIntensity();
      i++;
      j++;
    }
    else {
      if (peak_1->getPosition() < peak_2->getPosition()) {
        i++;
      }
      else {
        j++;
      }
    }
  }
  return dot_prod;
}

}  // namespace deconv_ms_util

}  // namespace toppic
