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

#ifndef TOPPIC_TOPFD_DECONV_DECONV_DATA_UTIL_HPP_
#define TOPPIC_TOPFD_DECONV_DECONV_DATA_UTIL_HPP_

#include "ms/spec/peak.hpp"
#include "topfd/deconv/deconv_data.hpp"

namespace toppic {

namespace deconv_data_util {

DeconvDataPtr getDataPtr(const PeakPtrVec &peak_list, double max_mass, 
                         int max_charge, double dp_window_size, 
                         bool estimate_min_inte, double sn_ratio);

};

}

#endif
