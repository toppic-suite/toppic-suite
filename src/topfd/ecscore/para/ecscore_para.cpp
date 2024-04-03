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

#include "common/util/file_util.hpp"
#include "common/base/mass_constant.hpp"
#include "topfd/ecscore/para/ecscore_para.hpp"

namespace toppic {

EcscorePara::EcscorePara(int frac_id, const std::string &file_name,
                         TopfdParaPtr para_ptr):
    frac_id_(frac_id),
    file_name_(file_name) {

  /// additional parameters
  para_max_charge_ = para_ptr->getMaxCharge();
  min_scan_num_ = para_ptr->getMinScanNum();
}

EcscorePara::EcscorePara(int frac_id, const std::string &file_name,
                         TopdiaParaPtr para_ptr, int ms_level):
        frac_id_(frac_id),
        file_name_(file_name) {

    /// additional parameters
    para_max_charge_ = para_ptr->getMaxCharge();
    if (ms_level == 1) {
      min_scan_num_ = para_ptr->getMs1MinScanNum();
      seed_env_inte_corr_tole_cutoff_ = para_ptr->getMs1SeedEnvInteCorrToleCutoff();
    }
    else {
      min_scan_num_ = para_ptr->getMs2MinScanNum();
      seed_env_inte_corr_tole_cutoff_ = para_ptr->getMs2SeedEnvInteCorrToleCutoff();
    }
}

} /* namespace */

