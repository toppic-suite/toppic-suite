//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#include "common/util/logger.hpp"
#include "ms/spec/raw_ms_util.hpp"
#include "topfd/spec/deconv_data_util.hpp"

namespace toppic {

namespace deconv_data_util {

DeconvDataPtr getDataPtr(const PeakPtrVec &peak_list, double max_mass, 
                         int max_charge, double window_size) {
  if (peak_list.size() == 0) return nullptr;
  double max_mz = raw_ms_util::findMaxPos(peak_list);
  if (max_mz > max_mass) {
    LOG_WARN("Max mz is too large: " << max_mz);
    return nullptr;
  }

  return std::make_shared<DeconvData>(peak_list, max_mass, 
                                      max_charge, window_size);
}


// generate deconvolution data using given max mass, max charge
DeconvDataPtr getDataPtr(const PeakPtrVec &peak_list, double spec_max_mass,
                         int spec_max_charge, double para_max_mass, 
                         int para_max_charge, double window_size) {
  if (spec_max_charge < 1) {
    LOG_INFO("Max charge < 1");
    spec_max_charge = para_max_charge;
  }
  if (spec_max_mass <= 0) {
    LOG_INFO("Max mass <= 0");
    spec_max_mass = para_max_mass;
  }
  if (spec_max_mass > para_max_mass) {
    LOG_WARN("Max mass is greater than default max mass " << spec_max_mass);
    spec_max_mass = para_max_mass;
  }
  double max_mz = raw_ms_util::findMaxPos(peak_list);
  if (max_mz > para_max_mass) {
    LOG_WARN("Max mz is too large: " << max_mz);
    return nullptr;
  }
  for (size_t i = 0; i < peak_list.size(); i++) {
    if (peak_list[i]->getPosition() < 0 || peak_list[i]->getIntensity() < 0) {
      LOG_WARN("mz intensity are negative values");
      return nullptr;
    }
  }
  return std::make_shared<DeconvData>(peak_list, spec_max_mass,
                                      spec_max_charge, window_size);
}

}

}  // namespace toppic
