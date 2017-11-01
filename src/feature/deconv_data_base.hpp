//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_DECONV_DATA_BASE_HPP_
#define PROT_FEATURE_DECONV_DATA_BASE_HPP_

#include <memory>
#include <vector>

#include "spec/peak.hpp"
#include "feature/feature_mng.hpp"
#include "feature/deconv_data.hpp"

namespace prot {

class DeconvDataBase {
 public:
  static DeconvDataPtr getDataPtr(PeakPtrVec &peak_list, FeatureMngPtr
                                  mng_ptr); 
  static DeconvDataPtr getDataPtr(PeakPtrVec &peak_list, double max_mass, int
                                  max_charge, FeatureMngPtr mng_ptr);
};

}

#endif
