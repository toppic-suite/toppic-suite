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

#include <numeric>
#include <algorithm>

#include "common/util/logger.hpp"

#include "ms/spec/peak_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/msmap/ms_map_row_header.hpp"

#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"
#include "topfd/ecscore/score/component_score.hpp"
#include "topfd/ecscore/score/comp_env_cnn_score.hpp"
#include "topfd/ecscore/score/onnx_ecscore.hpp"
#include "topfd/ecscore/score/ecscore.hpp"

namespace toppic {

ECScore::ECScore(EnvCollPtr env_coll_ptr, MsMapPtr matrix_ptr,
                 int score_id, double sn_ratio) {
  SeedEnvPtr seed_ptr = env_coll_ptr->getSeedPtr();
  MsMapRowHeaderPtrVec spec_list = matrix_ptr->getHeaderPtrList();

  EnvSetPtr env_set_ptr = env_coll_ptr->getSeedEnvSet();
  score_id_ = score_id;

  min_scan_ = env_coll_ptr->getStartSpecId();
  max_scan_ = env_coll_ptr->getEndSpecId();
  min_charge_ = env_coll_ptr->getMinCharge();
  max_charge_ = env_coll_ptr->getMaxCharge();
  mono_mass_ = seed_ptr->getMonoNeutralMass();
  rep_charge_ = seed_ptr->getCharge();
  rep_mz_ = seed_ptr->getMonoMz();

  abundance_ = env_coll_ptr->getIntensity();

  // convert seconds to minutes
  min_elution_time_ = spec_list[min_scan_]->getRt()/60;
  max_elution_time_ = spec_list[max_scan_]->getRt()/60; 
  int seed_spec_id = seed_ptr->getSpecId();
  apex_elution_time_ = spec_list[seed_spec_id]->getRt()/60;
  elution_length_ = max_elution_time_ - min_elution_time_; 

  EnvSetPtr seed_set_ptr = env_coll_ptr->getSeedEnvSet();
  double min_inte = matrix_ptr->getBaseInte() * sn_ratio;
  std::vector<std::vector<double>> theo_map 
    = seed_set_ptr->getScaledTheoIntes(min_inte);

  map_max_elution_time_ = spec_list[spec_list.size()-1]->getRt()/60;

  percent_matched_peaks_ = component_score::getMatchedPeakPercent(env_set_ptr, theo_map);
  intensity_correlation_ = component_score::getAggEnvCorr(env_set_ptr);
  top3_correlation_ = component_score::get3ScanCorr(env_set_ptr, seed_spec_id, min_scan_);
  even_odd_peak_ratio_ = component_score::getAggOddEvenPeakRatio(env_set_ptr);
  percent_consec_peaks_ = component_score::getConsecutivePeakPercent(env_set_ptr);
  num_theo_peaks_ = component_score::getTheoPeakNum(theo_map);
  mz_error_sum_ = component_score::getMzErrors(env_set_ptr);

  envcnn_score_ = comp_env_cnn_score::compEnvcnnScore(matrix_ptr, env_coll_ptr); 
  label_ = 0;

  std::vector<float> ecscore_input = getEcscoreInput(map_max_elution_time_);
  score_ = onnx_ecscore::predict(ecscore_input); 
}

std::vector<float> ECScore::getEcscoreInput(double max_elution_time) {
  std::vector<float> data;
  data.push_back(envcnn_score_); //1
  data.push_back(elution_length_ / max_elution_time * 2); //2
  data.push_back(percent_matched_peaks_); //3
  data.push_back(rep_charge_); //4
  data.push_back(top3_correlation_); //5
  data.push_back((max_charge_ - min_charge_) / 30.0); //6
  data.push_back(even_odd_peak_ratio_); //7
  return data;
}

}
