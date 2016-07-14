#include "base/logger.hpp"
#include "feature/deconv_data_base.hpp"
#include "feature/deconv_util.hpp"
#include "feature/deconv_one_sp.hpp"
#include "feature/match_env_util.hpp"
#include "feature/match_env_refine.hpp"
#include "feature/match_env_filter.hpp"
#include "feature/env_detect.hpp"
#include "feature/env_filter.hpp"
#include "feature/env_assign.hpp"
#include "feature/co_table.hpp"
#include "feature/dp_a.hpp"

namespace prot {

DeconvOneSp::DeconvOneSp(FeatureMngPtr mng_ptr): 
    mng_ptr_(mng_ptr) {
	}

void DeconvOneSp::setData(PeakPtrVec &peak_list) {
  data_ptr_ = DeconvDataBase::getDataPtr(peak_list, mng_ptr_);
}

void DeconvOneSp::setData(PeakPtrVec &peak_list, double max_mass, int max_charge) {
  data_ptr_ = DeconvDataBase::getDataPtr(peak_list, max_mass, max_charge, mng_ptr_);
}

void DeconvOneSp::run() {
  preprocess();
  LOG_DEBUG("preprocess complete");
  // envelope detection 
  MatchEnvPtr2D cand_envs = EnvDetect::getCandidate(data_ptr_, mng_ptr_);
  LOG_DEBUG("candidate complete");
  // envelope filter 
  EnvFilter::filter(cand_envs, data_ptr_, mng_ptr_);
  // assign envelopes to 1 Da windows 
  MatchEnvPtr2D win_envs = EnvAssign::assignWinEnv(cand_envs, data_ptr_,
                                                   mng_ptr_->env_num_per_window_);

  // prepare table for dp 
  if (mng_ptr_->check_double_increase_) {
    LOG_DEBUG("Generating coexisting table...");
    mng_ptr_->coexist_table_ = CoTable::initCoexistTable(win_envs,
                                                         mng_ptr_->score_error_tolerance_);
  }
  // dp 
  LOG_DEBUG("Generating Graph and DP...");
  DpA dp (data_ptr_, win_envs, mng_ptr_);
  MatchEnvPtrVec dp_envs = dp.getResult();

  result_envs_ = postprocess(dp_envs);
}

void DeconvOneSp::preprocess() {
  if (mng_ptr_->estimate_min_inte_) {
    PeakPtrVec peak_list = data_ptr_->getPeakList();
    std::vector<double> intes;
    for (size_t i = 0; i < peak_list.size(); i++) {
      intes.push_back(peak_list[i]->getIntensity());
    }
    double min_inte = DeconvUtil::getBaseLine(intes);
    mng_ptr_->setMinInte(min_inte);
  }
}

MatchEnvPtrVec DeconvOneSp::postprocess(MatchEnvPtrVec  &dp_envs) {
  // assign intensity
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  MatchEnvUtil::assignIntensity(peak_list, dp_envs);
  // refinement
  if (!mng_ptr_->output_multiple_mass_) {
    MatchEnvRefine::mzRefine(mng_ptr_, dp_envs);
  }

  MatchEnvPtrVec result_envs_ = dp_envs;
  if (mng_ptr_->do_final_filtering_) {
    result_envs_ = MatchEnvFilter::filter(dp_envs, data_ptr_->getMaxMass(),
                                          mng_ptr_);
  }

  if (mng_ptr_->keep_unused_peaks_) {
    MatchEnvUtil::addLowMassPeak(result_envs_, peak_list, mng_ptr_->mz_tolerance_);
  }
  // reassign intensity
  MatchEnvUtil::assignIntensity(peak_list, result_envs_);

  if (mng_ptr_->output_multiple_mass_) {
    // envelope detection 
    MatchEnvPtr2D cand_envs = EnvDetect::getCandidate(data_ptr_, mng_ptr_);
    // envelope filter 
    EnvFilter::multipleMassFilter(cand_envs, data_ptr_, mng_ptr_);
    result_envs_ = MatchEnvUtil::addMultipleMass(result_envs_, cand_envs,
                                                 mng_ptr_->multiple_min_mass_, 
                                                 mng_ptr_->multiple_min_charge_, 
                                                 mng_ptr_->multiple_min_ratio_);
  }

  return result_envs_;
}

}
