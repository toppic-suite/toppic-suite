//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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

#include <cmath>
#include <limits>

#include "common/util/logger.hpp"

#include "topfd/ecscore/env/ms_map_env_util.hpp"

namespace toppic {

namespace ms_map_env_util {

double pearsonr(std::vector<double> &X, std::vector<double> &Y) {
  int n = X.size();
  double sum_X = 0, sum_Y = 0, sum_XY = 0;
  double squareSum_X = 0, squareSum_Y = 0;
  for (int i = 0; i < n; i++) {
    sum_X = sum_X + X[i];
    sum_Y = sum_Y + Y[i];
    sum_XY = sum_XY + X[i] * Y[i];
    squareSum_X = squareSum_X + X[i] * X[i];
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
  }
  double corr = (double) (n * sum_XY - sum_X * sum_Y) /
    sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
  return corr;
}

double compTopThreeInteRatio(SeedEnvPtr seed_ptr, MsMapEnvPtr env_ptr) {
  double seed_inte = seed_ptr->compTopThreeInteSum();
  int ref_idx = seed_ptr->getReferIdx();
  double env_inte = env_ptr->compTopThreeInteSum(ref_idx);
  if (seed_inte == 0) {
    LOG_WARN("Empty peak list!");
    return 0;
  }
  return env_inte/seed_inte;
}

double compTopThreeInteRatio(SeedEnvPtr seed_ptr, std::vector<double> &inte_list) {
  double seed_inte = seed_ptr->compTopThreeInteSum();
  if (seed_inte == 0) {
    LOG_WARN("Empty peak list!");
    return 0;
  }
  size_t ref_idx = seed_ptr->getReferIdx();
  double env_inte = inte_list[ref_idx];
  if (ref_idx - 1 >= 0) {
    env_inte += inte_list[ref_idx-1];
  }
  if (ref_idx + 1 < inte_list.size()) {
    env_inte += inte_list[ref_idx+1];
  }
  return env_inte/seed_inte;
}

MsMapPeakPtr pickMsMapPeak(MsMapPtr ms_map_ptr, EnvPeakPtr seed_peak_ptr,
                           int sp_id, double mass_tol) {
  // get peaks within mass tolerance
  double max_inte = std::numeric_limits<double>::min();
  double mz = seed_peak_ptr->getPosition();
  int start_idx = ms_map_ptr->getColIndex(mz - mass_tol);
  if (start_idx < 0) {
      start_idx = 0;
  }
  int end_idx = ms_map_ptr->getColIndex(mz + mass_tol);
  if (end_idx >= ms_map_ptr->getColNum()) {
      end_idx = ms_map_ptr->getColNum() - 1;
  }
  MsMapPeakPtr result_peak = nullptr;
  for (int idx = start_idx; idx <= end_idx; idx++) {
    MsMapPeakPtrVec bin_peaks = ms_map_ptr->getBinPeakList(sp_id, idx);
    for (const auto& peak_ptr : bin_peaks) {
      double mass_diff = std::abs(mz - peak_ptr->getPosition());
      if ( mass_diff < mass_tol && peak_ptr->getIntensity() > max_inte) {
        result_peak = peak_ptr;
        max_inte = peak_ptr->getIntensity();
      }
    }
  }
  return result_peak;
}

MsMapEnvPtr getMatchMsMapEnv(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                             int sp_id, double mass_tol) {
  MsMapPeakPtrVec peak_list;
  EnvPeakPtrVec peaks = seed_ptr->getPeakPtrList();
  for (auto& seed_peak : peaks) {
    if (seed_peak != nullptr) {
      MsMapPeakPtr peak = pickMsMapPeak(ms_map_ptr, seed_peak, sp_id, mass_tol);
      peak_list.push_back(peak);
    }
    else {
      peak_list.push_back(nullptr);
    }
  }
  MsMapEnvPtr ms_map_env_ptr = std::make_shared<MsMapEnv>(sp_id, peak_list);
  return ms_map_env_ptr;
}

MsMapEnvPtr getMatchMsMapEnv(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                             int sp_id, double mass_tole,
                             double min_inte) {
  MsMapEnvPtr ms_map_env_ptr = getMatchMsMapEnv(ms_map_ptr, seed_ptr,
                                                sp_id, mass_tole);
  double inte_ratio = ms_map_env_util::compTopThreeInteRatio(seed_ptr, ms_map_env_ptr);
  ms_map_env_ptr->removeLowIntePeaks(seed_ptr, inte_ratio, min_inte);
  return ms_map_env_ptr;
}


}

}
