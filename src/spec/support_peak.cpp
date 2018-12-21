//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include "spec/support_peak.hpp"

namespace toppic {

SupportPeak::SupportPeak(DeconvPeakPtr deconv_peak_ptr, double offset,
                         double score, SPTypePtr peak_type_ptr): 
    deconv_peak_ptr_(deconv_peak_ptr),
    offset_(offset),
    score_(score),
    peak_type_ptr_(peak_type_ptr) {}

} /* namespace toppic */
