// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_COLL_ENV_COLL_UTIL_HPP_
#define TOPPIC_TOPFD_ECSCORE_ENV_COLL_ENV_COLL_UTIL_HPP_

#include <cstddef>
#include <vector>

#include "ms/feature/frac_feature.hpp"
#include "ms/msmap/ms_map.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "topfd/ecscore/env/seed_env.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"
#include "topfd/ecscore/para/ecscore_para.hpp"

namespace toppic {

namespace env_coll_util {

EnvCollPtr findEnvColl(const MsMapPtr& matrix_ptr, const SeedEnvPtr& seed_ptr,
                       const EcscoreParaPtr& para_ptr, double sn_ratio,
                       double split_ratio);

EnvCollPtr findEnvCollWithSingleEnv(const MsMapPtr& matrix_ptr,
                                    const SeedEnvPtr& seed_ptr,
                                    const EcscoreParaPtr& para_ptr,
                                    double sn_ratio, double split_ratio);

// The retention-time ranges of the features of an EnvCollPtrVec in
// contiguous arrays, for checkExistingFeatures, which scans the whole list
// for every candidate feature (88,000 candidates against up to 42,000
// features on 480_F3_02, 90 s spent chasing the feature and header objects;
// toppic_test experiment 36). sync appends the entries of the features
// added since the last call, refresh updates one feature after a merge,
// and clear must be called when the map is rebuilt (its retention times
// are copied here).
class EnvCollRtCache {
 public:
  void sync(const MsMapPtr& matrix_ptr, const EnvCollPtrVec& env_coll_list);
  void refresh(const MsMapPtr& matrix_ptr, const EnvCollPtr& coll_ptr,
               std::size_t idx);
  void clear() {
    start_rt_.clear();
    end_rt_.clear();
  }
  std::vector<double> start_rt_;
  std::vector<double> end_rt_;
};

// With a cache, the retention-time ranges of the listed features are read
// from it (same result, same first match) instead of from the feature and
// header objects.
bool checkExistingFeatures(const MsMapPtr& matrix_ptr,
                           const EnvCollPtr& coll_ptr,
                           EnvCollPtrVec& env_coll_list,
                           const EcscoreParaPtr& para_ptr, double sn_ratio,
                           EnvCollRtCache* cache = nullptr);

FracFeaturePtr getFracFeature(int feat_id, DeconvMsPtrVec& ms1_ptr_vec,
                              int frac_id, std::string& file_name,
                              const EnvCollPtr& coll_ptr,
                              const MsMapPtr& matrix_ptr, double snr);

}  // namespace env_coll_util
}  // namespace toppic

#endif
