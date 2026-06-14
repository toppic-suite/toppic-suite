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

#include "topfd/deconv/deconv_single_sp.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <mutex>
#include <ostream>

#include "common/util/logger.hpp"
#include "ms/env/env_detect.hpp"
#include "ms/env/env_filter.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/match_env_filter.hpp"
#include "ms/env/match_env_refine.hpp"
#include "ms/env/match_env_util.hpp"
#include "topfd/deconv/deconv_data_util.hpp"
#include "topfd/dp/co_table.hpp"
#include "topfd/dp/dp_a.hpp"
#include "topfd/dp/dp_assign.hpp"
#include "topfd/dp/dp_para.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

namespace {

// Serializes the debug dumps and numbers the spectra across worker threads.
std::mutex g_dp_env_mutex;
int g_dp_env_spec_count = 0;
int g_cand_env_spec_count = 0;

// Write one envelope's information: a summary line followed by its peaks.
void writeEnv(std::ostream& os, const MatchEnvPtr& match_env) {
  EnvPtr theo_env = match_env->getTheoEnvPtr();
  int peak_num = theo_env->getPeakNum();
  os << "  env mono_mass " << theo_env->getMonoNeutralMass() << " charge "
     << theo_env->getCharge() << " intensity " << theo_env->compInteSum()
     << " peak_num " << peak_num << "\n";
  for (int k = 0; k < peak_num; k++) {
    os << "    " << theo_env->getMz(k) << " " << theo_env->getInte(k) << "\n";
  }
}

// Dump the non-empty windowed candidate envelopes (win_envs) and the
// DP-selected envelopes (dp_envs) of one spectrum to win_envs.txt /
// dp_envs.txt. The first spectrum of the run truncates the files; later spectra
// append.
void outputDpEnvs(const MatchEnvPtr2D& win_envs, const MatchEnvPtrVec& dp_envs,
                  int ms_level) {
  std::lock_guard<std::mutex> lock(g_dp_env_mutex);
  std::ios::openmode mode =
      (g_dp_env_spec_count == 0) ? std::ios::out : std::ios::app;
  int spec_index = g_dp_env_spec_count++;

  std::ofstream win_out("win_envs.txt", mode);
  win_out << "# spectrum " << spec_index << " ms_level " << ms_level << "\n";
  for (size_t w = 0; w < win_envs.size(); w++) {
    if (win_envs[w].empty()) {
      continue;
    }
    win_out << "window " << w << " " << win_envs[w].size() << " envelopes\n";
    for (size_t i = 0; i < win_envs[w].size(); i++) {
      writeEnv(win_out, win_envs[w][i]);
    }
  }

  std::ofstream dp_out("dp_envs.txt", mode);
  dp_out << "# spectrum " << spec_index << " ms_level " << ms_level << " "
         << dp_envs.size() << " envelopes\n";
  for (size_t i = 0; i < dp_envs.size(); i++) {
    writeEnv(dp_out, dp_envs[i]);
  }
}

// Dump the non-null candidate envelopes (per peak and charge) of one spectrum
// to cand_envs.txt, before any filtering. The first spectrum truncates the
// file; later spectra append.
void outputCandEnvs(const MatchEnvPtr2D& cand_envs, int ms_level) {
  std::lock_guard<std::mutex> lock(g_dp_env_mutex);
  std::ios::openmode mode =
      (g_cand_env_spec_count == 0) ? std::ios::out : std::ios::app;
  int spec_index = g_cand_env_spec_count++;

  std::ofstream out("cand_envs.txt", mode);
  out << "# spectrum " << spec_index << " ms_level " << ms_level << "\n";
  for (size_t i = 0; i < cand_envs.size(); i++) {
    for (size_t j = 0; j < cand_envs[i].size(); j++) {
      if (cand_envs[i][j] != nullptr) {
        out << "peak " << i << " charge " << (j + 1) << "\n";
        writeEnv(out, cand_envs[i][j]);
      }
    }
  }
}

}  // namespace

DeconvSingleSp::DeconvSingleSp(const TopfdParaPtr& topfd_para_ptr,
                               PeakPtrVec& peak_list, int ms_level,
                               double max_mass, int max_charge) {
  topfd_para_ptr_ = topfd_para_ptr;
  env_para_ptr_ = std::make_shared<EnvPara>(topfd_para_ptr->getMzError());
  dp_para_ptr_ = std::make_shared<DpPara>(topfd_para_ptr->getMzError());
  ms_level_ = ms_level;
  double sn_ratio;
  if (ms_level == 1) {
    sn_ratio = topfd_para_ptr_->getMsOneSnRatio();
  } else {
    sn_ratio = topfd_para_ptr_->getMsTwoSnRatio();
  }
  data_ptr_ = deconv_data_util::getDataPtr(
      peak_list, max_mass, max_charge, dp_para_ptr_->dp_window_size_,
      topfd_para_ptr_->isEstimateMinInte(), sn_ratio);
}

void DeconvSingleSp::postprocess(MatchEnvPtrVec& dp_envs, int ms_level) {
  // assign intensity
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  match_env_util::assignIntensity(peak_list, dp_envs);

  // refinement
  match_env_refine::mzRefine(dp_envs);
  result_envs_ = dp_envs;

  if (topfd_para_ptr_->isKeepUnusedPeaks()) {
    // all added envelopes have charge 1, so we use charge 1 to get mz_tolerance
    double mz_tole = env_para_ptr_->getMzTolerance(1);
    result_envs_ =
        match_env_util::addUnusedMasses(result_envs_, peak_list, mz_tole);
  }

  // Obtain EnvCNN Scores for envelopes
  onnx_env_cnn::computeEnvScores(peak_list, result_envs_);
  // filtering by EnvCNN score
  if (ms_level == 2) {
    result_envs_ = match_env_filter::filterByEnvCnnScore(
        result_envs_, topfd_para_ptr_->getMs2EnvCnnScoreCutoff());
  }
  // filtering by estimated number of amino acids
  if (topfd_para_ptr_->isAANumBasedFilter()) {
    result_envs_ = match_env_filter::filterByAANum(
        result_envs_, data_ptr_->getMaxMass(),
        topfd_para_ptr_->isSortUseMsDeconv(), env_para_ptr_);
  }
  // sorting
  if (topfd_para_ptr_->isSortUseMsDeconv()) {
    std::sort(result_envs_.begin(), result_envs_.end(),
              MatchEnv::cmpMsdeconvScoreDec);
  } else {
    std::sort(result_envs_.begin(), result_envs_.end(),
              MatchEnv::cmpEnvcnnScoreDec);
  }
  // reassign intensity
  match_env_util::assignIntensity(peak_list, result_envs_);

  // if output multiple masses sharing the same envelope
  if (topfd_para_ptr_->isOutputMultipleMass()) {
    MatchEnvPtr2D cand_envs = env_detect::getCandidateEnv(
        peak_list, data_ptr_->getMaxCharge(), data_ptr_->getMaxMass(),
        data_ptr_->getMinInte(), data_ptr_->getMinRefInte(), env_para_ptr_);
    // envelope filter
    env_filter::multipleMassFilter(cand_envs, env_para_ptr_);
    result_envs_ =
        match_env_util::addMultipleMass(result_envs_, cand_envs, env_para_ptr_);
  }
}

MatchEnvPtrVec DeconvSingleSp::deconv() {
  if (data_ptr_ == nullptr) {
    return result_envs_;
  }
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  // envelope detection
  MatchEnvPtr2D cand_envs = env_detect::getCandidateEnv(
      peak_list, data_ptr_->getMaxCharge(), data_ptr_->getMaxMass(),
      data_ptr_->getMinInte(), data_ptr_->getMinRefInte(), env_para_ptr_);
  LOG_DEBUG("candidate complete");

  if (topfd_para_ptr_->isOutputDpEnvs()) {
    outputCandEnvs(cand_envs, ms_level_);
  }

  // envelope filter
  env_filter::filter(cand_envs, peak_list, env_para_ptr_);
  // prepare for dp
  MatchEnvPtr2D win_envs = dp_assign::assignWinEnv(
      cand_envs, data_ptr_->getWinNum(), data_ptr_->getWinIdVec(),
      dp_para_ptr_->env_num_per_win_);

  // prepare table for dp
  if (dp_para_ptr_->check_double_increase_) {
    LOG_DEBUG("Generating coexisting table...");
    dp_para_ptr_->coexist_table_ =
        CoTable::initCoexistTable(win_envs, dp_para_ptr_->mz_tolerance_);
  }
  // dp
  LOG_DEBUG("Generating Graph and DP...");
  DpA dp(data_ptr_, win_envs, dp_para_ptr_, dp_para_ptr_->mz_tolerance_);
  MatchEnvPtrVec dp_envs = dp.getResult();

  if (topfd_para_ptr_->isOutputDpEnvs()) {
    outputDpEnvs(win_envs, dp_envs, ms_level_);
  }

  postprocess(dp_envs, ms_level_);

  return result_envs_;
}

}  // namespace toppic
