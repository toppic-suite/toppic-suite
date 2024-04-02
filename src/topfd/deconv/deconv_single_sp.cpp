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

#include <algorithm>

#include "common/util/logger.hpp"

#include "ms/env/env_para.hpp"
#include "ms/env/env_detect.hpp"
#include "ms/env/env_filter.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/match_env_refine.hpp"
#include "ms/env/match_env_filter.hpp"

#include "topfd/dp/dp_para.hpp"
#include "topfd/dp/dp_assign.hpp"
#include "topfd/dp/co_table.hpp"
#include "topfd/dp/dp_a.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp" 
#include "topfd/deconv/deconv_data_util.hpp"
#include "topfd/deconv/deconv_single_sp.hpp"

namespace toppic {

DeconvSingleSp::DeconvSingleSp(TopfdParaPtr topfd_para_ptr, PeakPtrVec &peak_list,
                               int ms_level, double max_mass, int max_charge) {
  topfd_para_ptr_ = topfd_para_ptr; 
  env_para_ptr_ = std::make_shared<EnvPara>(topfd_para_ptr->getMzError());
  dp_para_ptr_ = std::make_shared<DpPara>(topfd_para_ptr->getMzError());
  double sn_ratio;
  if (ms_level == 1) {
    sn_ratio = topfd_para_ptr_->getMsOneSnRatio();
  }
  else {
    sn_ratio = topfd_para_ptr_->getMsTwoSnRatio();
  }
  data_ptr_ = deconv_data_util::getDataPtr(peak_list, max_mass, max_charge, 
                                           dp_para_ptr_->dp_window_size_,  
                                           topfd_para_ptr_->isEstimateMinInte(),
                                           sn_ratio);
}

void DeconvSingleSp::postprocess(MatchEnvPtrVec &dp_envs) {
  // assign intensity
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  match_env_util::assignIntensity(peak_list, dp_envs);
  // refinement
  if (!topfd_para_ptr_->isOutputMultipleMass()) {
    match_env_refine::mzRefine(dp_envs);
  }

  // Obtain EnvCNN Scores for envelopes
  onnx_env_cnn::computeEnvScores(peak_list, dp_envs); 
  if (topfd_para_ptr_->isUseMsDeconv()) {
    std::sort(dp_envs.begin(), dp_envs.end(), MatchEnv::cmpMsdeconvScoreDec);
  }
  else {
    std::sort(dp_envs.begin(), dp_envs.end(), MatchEnv::cmpEnvcnnScoreDec);
  }

  // filtering
  if (topfd_para_ptr_->isDoFinalFiltering()) {
    result_envs_ = match_env_filter::filter(dp_envs, data_ptr_->getMaxMass(), 
                                            topfd_para_ptr_->isUseMsDeconv(), 
                                            env_para_ptr_);
  }
  else {
    result_envs_ = dp_envs;
  }
  if (topfd_para_ptr_->isKeepUnusedPeaks()) {
    // all added envelopes have charge 1, so we use charge 1 to get mz_tolerance
//    double mz_tole = env_para_ptr_->getMzTolerance(1);
//    match_env_util::addUnusedMasses(result_envs_, peak_list, mz_tole);
    match_env_util::addUnusedMasses(result_envs_, peak_list,
                                    env_para_ptr_->getMzTolerance());
  }
  // reassign intensity
  match_env_util::assignIntensity(peak_list, result_envs_);

  // if output multiple masses sharing the same envelope
  if (topfd_para_ptr_->isOutputMultipleMass()) {
    MatchEnvPtr2D cand_envs = env_detect::getCandidateEnv(peak_list, 
                                                          data_ptr_->getMaxCharge(), 
                                                          data_ptr_->getMaxMass(), 
                                                          data_ptr_->getMinInte(), 
                                                          data_ptr_->getMinRefInte(),
                                                          env_para_ptr_);
    // envelope filter
    env_filter::multipleMassFilter(cand_envs, env_para_ptr_);
    result_envs_ = match_env_util::addMultipleMass(result_envs_, cand_envs, env_para_ptr_);
  }
}

MatchEnvPtrVec DeconvSingleSp::deconv() {
  if (data_ptr_ == nullptr) {
    return result_envs_;
  }
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  // envelope detection
  MatchEnvPtr2D cand_envs = env_detect::getCandidateEnv(peak_list,   
                                                        data_ptr_->getMaxCharge(),
                                                        data_ptr_->getMaxMass(),
                                                        data_ptr_->getMinInte(), 
                                                        data_ptr_->getMinRefInte(),
                                                        env_para_ptr_);
  LOG_DEBUG("candidate complete");
  // envelope filter
  env_filter::filter(cand_envs, peak_list, env_para_ptr_);
  // prepare for dp
  MatchEnvPtr2D win_envs = dp_assign::assignWinEnv(cand_envs, data_ptr_->getWinNum(),
                                                   data_ptr_->getWinIdVec(),
                                                   dp_para_ptr_->env_num_per_win_);

  // prepare table for dp
  if (dp_para_ptr_->check_double_increase_) {
    LOG_DEBUG("Generating coexisting table...");
    dp_para_ptr_->coexist_table_ = CoTable::initCoexistTable(win_envs,
                                                             dp_para_ptr_->mz_tolerance_);
  }
  // dp
  LOG_DEBUG("Generating Graph and DP...");
  DpA dp(data_ptr_, win_envs, dp_para_ptr_, dp_para_ptr_->mz_tolerance_); 
  MatchEnvPtrVec dp_envs = dp.getResult();

  postprocess(dp_envs);
  return result_envs_;
}

}  // namespace toppic

