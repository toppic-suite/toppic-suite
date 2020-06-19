//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <src/envcnn/score.hpp>
#include "common/util/logger.hpp"
#include "ms/spec/baseline_util.hpp"
#include "topfd/spec/deconv_data_util.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/match_env_refine.hpp"
#include "ms/env/match_env_filter.hpp"
#include "ms/env/env_detect.hpp"
#include "ms/env/env_filter.hpp"
#include "ms/env/env_assign.hpp"
#include "topfd/dp/co_table.hpp"
#include "topfd/dp/dp_a.hpp"
#include "topfd/deconv/deconv_one_sp.hpp"

//#include "ms/env/env_rescore.hpp"

namespace toppic {

void DeconvOneSp::setData(PeakPtrVec &peak_list) {
  data_ptr_ = deconv_data_util::getDataPtr(peak_list, env_para_ptr_->max_mass_,
                                           env_para_ptr_->max_charge_, env_para_ptr_->window_size_);
}

void DeconvOneSp::setData(PeakPtrVec &peak_list, double spec_max_mass, int spec_max_charge) {
  data_ptr_ = deconv_data_util::getDataPtr(peak_list, spec_max_mass, spec_max_charge, 
                                           env_para_ptr_->max_mass_, env_para_ptr_->max_charge_,
                                           env_para_ptr_->window_size_);
}

void DeconvOneSp::run() {
  if (data_ptr_ == nullptr) {
    result_envs_.clear();
    return;
  }

  preprocess();
  LOG_DEBUG("preprocess complete");
  // envelope detection
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  MatchEnvPtr2D cand_envs = EnvDetect::getCandidate(peak_list, 
                                                    data_ptr_->getMaxCharge(),
                                                    data_ptr_->getMaxMass(),
                                                    env_para_ptr_);
  LOG_DEBUG("candidate complete");
  // envelope filter
  EnvFilter::filter(cand_envs, peak_list, env_para_ptr_);
  // envelope rescoring
  //if (ms_level_ > 1) {
    //EnvRescore::rescore(cand_envs, env_para_ptr_->env_rescore_para_);
  //}
  // assign envelopes to 1 Da windows
  MatchEnvPtr2D win_envs = env_assign::assignWinEnv(cand_envs, data_ptr_->getWinNum(),
                                                   data_ptr_->getWinIdVec(),
                                                   env_para_ptr_->env_num_per_window_);

  // prepare table for dp
  if (dp_para_ptr_->check_double_increase_) {
    LOG_DEBUG("Generating coexisting table...");
    dp_para_ptr_->coexist_table_ = CoTable::initCoexistTable(win_envs,
                                                             env_para_ptr_->getScoreErrorTolerance());
  }
  // dp
  LOG_DEBUG("Generating Graph and DP...");
  DpA dp(data_ptr_, win_envs, dp_para_ptr_, env_para_ptr_->score_error_tolerance_);
  MatchEnvPtrVec dp_envs = dp.getResult();

  result_envs_ = postprocess(dp_envs);
}

void DeconvOneSp::preprocess() {
  if (env_para_ptr_->estimate_min_inte_) {
    PeakPtrVec peak_list = data_ptr_->getPeakList();
    std::vector<double> intes;
    for (size_t i = 0; i < peak_list.size(); i++) {
      intes.push_back(peak_list[i]->getIntensity());
    }
    double min_inte = baseline_util::getBaseLine(intes);
    env_para_ptr_->setMinInte(min_inte, ms_level_);
  }
}

MatchEnvPtrVec DeconvOneSp::postprocess(MatchEnvPtrVec  &dp_envs) {
  // assign intensity
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  match_env_util::assignIntensity(peak_list, dp_envs);
  // refinement
  if (!env_para_ptr_->output_multiple_mass_) {
    match_env_refine::mzRefine(dp_envs);
  }

    /////////////////////////////////////// EnvCNN Changes ////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /* Obtain EnvCNN Prediction Score for MS/MS envelopes */

//    if (ms_level_ > 1){
//        result_envs_ = MatchEnvFilterCNN::filter_using_cnn(dp_envs, peak_list, model_);
//    }

    ///////////////////////////////////////////////////////////////////////////////////////////////

    // filtering
  if (env_para_ptr_->do_final_filtering_) {
    result_envs_ = MatchEnvFilter::filter(dp_envs, data_ptr_->getMaxMass(),
                                          env_para_ptr_);
  }
  else {
    result_envs_ = dp_envs;
  }
  if (env_para_ptr_->keep_unused_peaks_) {
    match_env_util::addLowMassPeak(result_envs_, peak_list, env_para_ptr_->getMzTolerance());
  }
  // reassign intensity
  match_env_util::assignIntensity(peak_list, result_envs_);

  if (env_para_ptr_->output_multiple_mass_) {
    // envelope detection
    PeakPtrVec peak_list = data_ptr_->getPeakList();
    MatchEnvPtr2D cand_envs = EnvDetect::getCandidate(peak_list, data_ptr_->getMaxCharge(), 
                                                      data_ptr_->getMaxMass(), env_para_ptr_);
    // envelope filter
    EnvFilter::multipleMassFilter(cand_envs, env_para_ptr_);
    result_envs_ = match_env_util::addMultipleMass(result_envs_, cand_envs,
                                                 env_para_ptr_->multiple_min_mass_,
                                                 env_para_ptr_->multiple_min_charge_,
                                                 env_para_ptr_->multiple_min_ratio_);
  }

  return result_envs_;
}

}  // namespace toppic
