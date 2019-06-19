//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "spec/peak.hpp"
#include "spec/raw_ms_util.hpp"
#include "feature/peak_cluster.hpp"

namespace toppic {

int PeakCluster::even_charge_idx_ = 0;
int PeakCluster::odd_charge_idx_ = 1;
double PeakCluster::win_size_ = 6.0;

double getBcDistance(std::vector<double> &v1, std::vector<double> &v2) {
  if (v1.size() != v2.size() || v1.size() == 0) {
    LOG_ERROR("Two vectors have different sizes!");
    exit(EXIT_FAILURE);
  }

  double s1 = 0.0;
  double s2 = 0.0;
  for (size_t i = 0; i < v1.size(); i++) {
    s1 += v1[i];
    s2 += v2[i];
  }

  if (!(s1 > 0) || !(s2 > 0)) return 10.0;

  double bc = 0.0;
  for (size_t i = 0; i < v1.size(); i++) {
    double p = v1[i]/ s1;
    double q = v2[i]/ s2;
    bc += std::sqrt(p * q); 
  }
  if (bc > std::exp(-10)) {
    return -std::log(bc);
  }
  else {
    return 10;
  }
}


double getPearsonCorr(std::vector<double> &v1, std::vector<double> &v2) {
  if (v1.size() != v2.size() || v1.size() == 0) {
    LOG_ERROR("Two vectors have different sizes or empty vector!");
    exit(EXIT_FAILURE);
  }
  if (v1.size() == 1) {
    return 1.0;
  }
  // Compute means
  double m1 = 0.0;
  double m2 = 0.0;
  for (size_t i = 0; i < v1.size(); i++) {
    m1 += v1[i];
    m2 += v2[i];
  }
  m1 /= v1.size();
  m2 /= v2.size();

  // compute Pearson correlation
  double cov = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;

  for (size_t i = 0; i < v1.size(); i++) {
    double d1 = v1[i] - m1;
    double d2 = v2[i] - m2;
    cov += d1 * d2;
    s1 += d1 * d1;
    s2 += d2 * d2;
  }

  if (s1 <= 0 || s2 <= 0) return 0;

  return cov < 0 ? 0 : cov / std::sqrt(s1 * s2);
}

PeakCluster::PeakCluster(MatchEnvPtr match_env) {
  theo_env_ = match_env->getTheoEnvPtr();
  rep_mass_ = theo_env_->getMonoNeutralMass();
  rep_charge_ = theo_env_->getCharge();

  int peak_num = theo_env_->getPeakNum();
  rep_summed_intensities_.resize(peak_num, 0.0);

  clearScores();

  flag_ = 0;
  init_score_ = false;
  smoother_ = std::make_shared<SavitzkyGolay>(9, 2);
}

void PeakCluster::addEnvelopes(int min_charge, int max_charge, 
                               int min_ms1_id, int max_ms1_id, 
                               int scan_begin, int scan_end,
                               RealEnvPtrVec envs) {
  int row_num = max_charge - min_charge + 1;
  int col_num = max_ms1_id - min_ms1_id + 1;
  
  min_charge_ = min_charge;
  max_charge_ = max_charge;

  min_ms1_id_ = min_ms1_id;
  max_ms1_id_ = max_ms1_id;

  scan_begin_ = scan_begin;
  scan_end_ = scan_end;

  real_envs_.resize(row_num);

  for (int i = 0; i < row_num; i++) {
    real_envs_[i].resize(col_num);
  }

  for (size_t i = 0; i < envs.size(); i++) {
    int row = envs[i]->getCharge() - min_charge_;
    int col = envs[i]->getSpId() - min_ms1_id;
    if (row >= 0 && row < row_num && col >= 0 && col < col_num) {
      real_envs_[row][col] = envs[i];
    }
  }
}

void PeakCluster::clearScores() {
  inte_distr_.resize(2, 0.0);
  best_corr_scores_.resize(2, 0.0);
  best_inte_scores_.resize(2, 0.0);
  best_dist_scores_.resize(2, 1.0);

  best_charges_.resize(2, 0);
  env_dist_scores_.resize(2, 1.0);
  env_corr_scores_.resize(2, 0.0);
  env_inte_scores_.resize(2, 0.0);
  xic_corr_between_best_charges_.resize(2, 0.0);
}

void PeakCluster::updateScore(RawMsPtrVec spec_list, bool check_pvalue) {
  int row_num = max_charge_ - min_charge_ + 1;
  int col_num = max_ms1_id_ - min_ms1_id_ + 1;
  int ref_idx = theo_env_->getReferIdx(); 

  clearScores();

  std::vector<double> best_charge_dists{10.0, 10.0};

  // sum up peak intensities
  int peak_num = theo_env_->getPeakNum();
  std::vector<double> theo_intensities = theo_env_->getIntensities();
  std::vector<double> summed_intensities(peak_num, 0);

  int xic_len = col_num + 18;
  int xic_start_idx = 9;

  std::vector<std::vector<double>> xic2(2);
  xic2[0].resize(xic_len, 0.0);
  xic2[1].resize(xic_len, 0.0);

  std::vector<std::vector<double>> charge_xic(row_num);

  double tmp_best_bc_dist = 10.0;
  double rep_env_bc_dist = 10.0;
  RealEnvPtr rep_env(nullptr);

  double rep_env_bc_dist_2 = 10.0;
  RealEnvPtr rep_env_2(nullptr);

  std::vector<double> tmp_best_dist_scores{10.0, 10.0};
  std::vector<double> tmp_best_inte_scores(2, 0.0);
  std::vector<double> tmp_best_corr_scores(2, 0.0);

  for (int i = 0; i < row_num; i++) {
    int charge = i + min_charge_;
    double ref_neutral_mass = theo_env_->getRefNeutralMass();
    double ref_mz = Peak::compMz(ref_neutral_mass, charge); 
    std::fill(summed_intensities.begin(), summed_intensities.end(), 0.0);

    std::vector<double> cur_xic(xic_len, 0.0);

    int charge_idx = (charge % 2 == 0) ? even_charge_idx_:odd_charge_idx_;
    // summed_most_abu_isotope_intensity
    double summed_iso_high_inte = 0.0;
    //summed_referenc_intensity
    double summed_win_high_inte = 0.0;

    for (int j = 0; j < col_num; j++) {
      RealEnvPtr env = real_envs_[i][j];
      if (env == nullptr) continue;
      
      // sum peak inte
      for (int k = 0; k < peak_num; k++) {
        summed_intensities[i] += env->getIntensity(k);
      }

      int ms1_id = min_ms1_id_ + j;
      PeakPtrVec all_peaks = spec_list[ms1_id]->getPeakPtrVec();
      PeakPtrVec win_peaks = raw_ms_util::getPeaksInWindow(all_peaks, ref_mz, win_size_);
      double win_high_inte = raw_ms_util::getHighestPeakInte(win_peaks);
      double win_median_inte = raw_ms_util::getMedianPeakInte(win_peaks);
       
      if (env->isExist(ref_idx)) {
        summed_iso_high_inte += env->getIntensity(ref_idx);
        summed_win_high_inte += win_high_inte;
      }
      double env_inte_sum = env->getIntensitySum();
      inte_distr_[charge_idx] += env_inte_sum; 

      std::vector<double> real_intensities = env->getIntensities();

      double new_bc_dist = getBcDistance(theo_intensities, real_intensities);
      double new_corr = getPearsonCorr(theo_intensities, real_intensities);

      bool good_env = (new_bc_dist < 0.07 || new_corr > 0.7);
      if (good_env) {
        xic2[charge_idx][xic_start_idx + j] += env_inte_sum;
        charge_xic[i][xic_start_idx+j] += env_inte_sum;
      }

      bool level_one_env = true;
      bool level_two_env = true;
      /*
      if (pValueCheck)
      {
        var poissonPValue = localWin.GetPoissonTestPValue(envelope.Peaks, TheoreticalEnvelope.Size);
        var rankSumPValue = localWin.GetRankSumTestPValue(envelope.Peaks, TheoreticalEnvelope.Size);
        levelOneEnvelope = (rankSumPValue < 0.01 && poissonPValue < 0.01);
        //levelTwoEnvelope = (rankSumPValue < 0.05 || poissonPValue < 0.05);
      }
      */
      if (level_one_env ) {
        if (new_bc_dist < best_dist_scores_[charge_idx]) {
          best_dist_scores_[charge_idx] = new_bc_dist;
          double new_inte_score = 1.0;
          if (win_median_inte > 0.0) {
            new_inte_score = env->getReferIntensity()/win_high_inte; 
          }
          best_inte_scores_[charge_idx] = std::max(best_inte_scores_[charge_idx], new_inte_score);
        }
        best_corr_scores_[charge_idx] = std::max(best_corr_scores_[charge_idx], new_corr);

        if (new_bc_dist < rep_env_bc_dist) {
          rep_env_bc_dist = new_bc_dist;
          rep_env = env;
        }
        /*
        // in the initial scoring, classify major and minor envelopes
        if (!_initScore && goodEnvelope) envelope.GoodEnough = true;
        */
      }
      if (level_two_env) {
        if (new_bc_dist < tmp_best_dist_scores[charge_idx]) {
          tmp_best_dist_scores[charge_idx] = new_bc_dist;
          double new_inte_score = 1.0;
          if (win_median_inte > 0.0) {
            new_inte_score = env->getReferIntensity()/win_high_inte; 
          }
          tmp_best_inte_scores[charge_idx] = std::max(tmp_best_inte_scores[charge_idx], new_inte_score);
        }
        tmp_best_corr_scores[charge_idx] = std::max(tmp_best_corr_scores[charge_idx], new_corr);

        if (new_bc_dist < rep_env_bc_dist_2) {
          rep_env_bc_dist_2 = new_bc_dist;
          rep_env_2 = env;
        }
      }

      double bc_dist = getBcDistance(theo_intensities, summed_intensities);
      env_dist_scores_[charge_idx] = std::min(env_dist_scores_[charge_idx], bc_dist);
      double pc = getPearsonCorr(theo_intensities, summed_intensities);
      env_corr_scores_[charge_idx] = std::max(env_corr_scores_[charge_idx], pc);

      if (best_charges_[charge_idx] < 1 || bc_dist < best_charge_dists[charge_idx]) {
        best_charges_[charge_idx] = charge;
        best_charge_dists[charge_idx] = bc_dist;
        if (summed_win_high_inte > 0.0) {
          env_inte_scores_[charge_idx] = summed_iso_high_inte/summed_win_high_inte;
        }
      }

      if (bc_dist < tmp_best_bc_dist) {
        tmp_best_bc_dist = bc_dist;
        rep_summed_intensities_ = summed_intensities;
      }
    }
  }

  // when good envelope is observed at only even charge...
  if (best_corr_scores_[0] > 0.7 && best_corr_scores_[1] < 0.5) {
    int idx = 1;
    best_corr_scores_[idx] = tmp_best_corr_scores[idx];
    best_inte_scores_[idx] = tmp_best_inte_scores[idx];
    best_dist_scores_[idx] = tmp_best_dist_scores[idx];
  }

  // when good envelope is observed at only odd charge...
  if (best_corr_scores_[1] > 0.7 && best_corr_scores_[0] < 0.5) {
    int idx = 0;
    best_corr_scores_[idx] = tmp_best_corr_scores[idx];
    best_inte_scores_[idx] = tmp_best_inte_scores[idx];
    best_dist_scores_[idx] = tmp_best_dist_scores[idx];
  }

  // normalize intensities
  double s = inte_distr_[0] + inte_distr_[1];
  if (s > 0) {
    inte_distr_[0] = inte_distr_[0] / s;
    inte_distr_[1] = inte_distr_[1] / s;
  }

  if (col_num > 1) {
    int even_best_charge = best_charges_[even_charge_idx_] - min_charge_;
    int odd_best_charge = best_charges_[odd_charge_idx_] - min_charge_;
    std::vector<double> v1 = smoother_->smooth(charge_xic[even_best_charge]);
    std::vector<double> v2 = smoother_->smooth(charge_xic[odd_best_charge]);
    xic_corr_between_best_charges_[0] = getPearsonCorr(v1, v2); 
    v1 = smoother_->smooth(xic2[even_charge_idx_]);
    v2 = smoother_->smooth(xic2[odd_charge_idx_]);
    xic_corr_between_best_charges_[1] = getPearsonCorr(v1, v2);
  }

  if (rep_env == nullptr && rep_env_2 != nullptr) {
    rep_env = rep_env_2;
  }
  if (rep_env != nullptr) {
    rep_charge_ = rep_env->getCharge();
    rep_mz_ = rep_env->getMonoMz();
    rep_ms1_id_ = rep_env->getSpId();
  }

  init_score_ = true;
}

}

