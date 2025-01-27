//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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
#include "ms/spec/peak_list_util.hpp"
#include "topfd/deconv/deconv_data_util.hpp"

namespace toppic {

namespace deconv_data_util {

DeconvDataPtr getDataPtr(const PeakPtrVec &peak_list, 
                         double max_mass, int max_charge, 
                         double dp_window_size,
                         bool estimate_min_inte, 
                         double sn_ratio) {
  if (peak_list.size() == 0) return nullptr;
  double max_mz = peak_list_util::findMaxPos(peak_list);
  if (max_mz > max_mass) {
    LOG_INFO("Max mz is too large: " << max_mz);
  }

  return std::make_shared<DeconvData>(peak_list, max_mass, 
                                      max_charge, dp_window_size, 
                                      estimate_min_inte, sn_ratio);
}


}

}  // namespace toppic
