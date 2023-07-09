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

#include <algorithm>

#include "common/util/logger.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/match_env_refine.hpp"
#include "ms/env/match_env_filter.hpp"
#include "ms/env/env_detect.hpp"
#include "ms/env/env_filter.hpp"
#include "ms/env/env_assign.hpp"
#include "topfd/spec/deconv_data_util.hpp"
#include "topfd/dp/co_table.hpp"
#include "topfd/dp/dp_a.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp" 
#include "topfd/deconv/deconv_one_sp.hpp"

namespace toppic {

void DeconvOneSp::setData(PeakPtrVec &peak_list) {
  data_ptr_ = deconv_data_util::getDataPtr(peak_list, topfd_para_ptr_->getMaxMass(),
                                           topfd_para_ptr_->getMaxCharge(), topfd_para_ptr_->getPrecWindow());
}

void DeconvOneSp::setData(PeakPtrVec &peak_list, double spec_max_mass, int spec_max_charge) {
  data_ptr_ = deconv_data_util::getDataPtr(peak_list, spec_max_mass, spec_max_charge, 
                                           topfd_para_ptr_->getMaxMass(), topfd_para_ptr_->getMaxCharge(),
                                           topfd_para_ptr_->getDpWindowSize());
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
  MatchEnvPtr2D cand_envs = env_detect::getCandidateEnv(peak_list, 
                                                        data_ptr_->getMaxCharge(),
                                                        data_ptr_->getMaxMass(),
                                                        data_ptr_->getMinInte(), 
                                                        data_ptr_->getMinRefInte(),
                                                        env_para_ptr_);
  LOG_DEBUG("candidate complete");
  // envelope filter
  env_filter::filter(cand_envs, peak_list, env_para_ptr_);
  // envelope rescoring
  //if (ms_level_ > 1) {
    //EnvRescore::rescore(cand_envs, env_para_ptr_->env_rescore_para_);
  //}
  // assign envelopes to 1 Da windows
  MatchEnvPtr2D win_envs = env_assign::assignWinEnv(cand_envs, data_ptr_->getWinNum(),
                                                    data_ptr_->getWinIdVec(),
                                                    topfd_para_ptr_->getEnvNumPerWin());

  // prepare table for dp
  if (dp_para_ptr_->check_double_increase_) {
    LOG_DEBUG("Generating coexisting table...");
    dp_para_ptr_->coexist_table_ = CoTable::initCoexistTable(win_envs,
                                                             env_para_ptr_->getScoreErrorTolerance());
  }
  // dp
  LOG_DEBUG("Generating Graph and DP...");
  DpA dp(data_ptr_, win_envs, dp_para_ptr_, env_para_ptr_->getScoreErrorTolerance());
  MatchEnvPtrVec dp_envs = dp.getResult();

  result_envs_ = postprocess(dp_envs);
}

void DeconvOneSp::preprocess() {
  double min_inte = 0;
  double min_ref_inte = 0;
  if (topfd_para_ptr_->isEstimateMinInte()) {
    PeakPtrVec peak_list = data_ptr_->getPeakList();
    std::vector<double> intes;
    for (size_t i = 0; i < peak_list.size(); i++) {
      intes.push_back(peak_list[i]->getIntensity());
    }
    min_inte = baseline_util::getBaseLine(intes);
    if (ms_level_ == 1) {
      min_ref_inte = min_inte * topfd_para_ptr_->getMsOneSnRatio();
    }
    else {
      min_ref_inte = min_inte * topfd_para_ptr_->getMsTwoSnRatio();
    }
  }
  data_ptr_->setMinInte(min_inte); 
  data_ptr_->setMinRefInte(min_ref_inte);
}

MatchEnvPtrVec DeconvOneSp::postprocess(MatchEnvPtrVec &dp_envs) {
  // assign intensity
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  match_env_util::assignIntensity(peak_list, dp_envs);
  // refinement
  if (!topfd_para_ptr_->isOutputMultipleMass()) {
    match_env_refine::mzRefine(dp_envs);
  }

  // Obtain EnvCNN Prediction Score for MS/MS envelopes
  if (topfd_para_ptr_->isUseEnvCnn() && ms_level_ != 1){
    //env_cnn::computeEnvScores(dp_envs, peak_list);
    onnx_env_cnn::computeEnvScores(peak_list, dp_envs); 
    std::sort(dp_envs.begin(), dp_envs.end(), MatchEnv::cmpEnvcnnScoreDec);
  }

  // filtering
  if (topfd_para_ptr_->isDoFinalFiltering()) {
    result_envs_ = match_env_filter::filter(dp_envs, data_ptr_->getMaxMass(),
                                            topfd_para_ptr_->isUseEnvCnn(), env_para_ptr_);
  }
  else {
    result_envs_ = dp_envs;
  }
  if (topfd_para_ptr_->isKeepUnusedPeaks()) {
    match_env_util::addUnusedMasses(result_envs_, peak_list, env_para_ptr_->getMzTolerance());
  }
  // reassign intensity
  match_env_util::assignIntensity(peak_list, result_envs_);

  if (topfd_para_ptr_->isOutputMultipleMass()) {
    // envelope detection
    PeakPtrVec peak_list = data_ptr_->getPeakList();
    MatchEnvPtr2D cand_envs = env_detect::getCandidateEnv(peak_list, data_ptr_->getMaxCharge(), 
                                                          data_ptr_->getMaxMass(), 
                                                          data_ptr_->getMinInte(), 
                                                          data_ptr_->getMinRefInte(),
                                                          env_para_ptr_);
    // envelope filter
    env_filter::multipleMassFilter(cand_envs, env_para_ptr_);
    result_envs_ = match_env_util::addMultipleMass(result_envs_, cand_envs, env_para_ptr_);
  }

  return result_envs_;
}

}  // namespace toppic
