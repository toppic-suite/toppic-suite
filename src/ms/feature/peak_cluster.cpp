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

#include <cmath>

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>

#include "common/util/logger.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/env_peak.hpp"
#include "ms/spec/raw_ms_util.hpp"
#include "ms/feature/peak_cluster.hpp"

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

double getPoissonPValue(PeakPtrVec &win_peaks, double win_size, 
                        std::vector<double> &env_peak_intensities) {
  double win_high_inte = raw_ms_util::getHighestPeakInte(win_peaks);
  double inte_thresh = win_high_inte  * 0.1;
  int match_peak_num = 0;
  for (size_t i = 0; i < env_peak_intensities.size(); i++) {
    if (env_peak_intensities[i] > inte_thresh) {
      match_peak_num++;
    }
  }
  int inte_peak_num = 0;
  for (size_t i = 0; i < win_peaks.size(); i++) {
    if (win_peaks[i]->getIntensity() > inte_thresh) {
      inte_peak_num++;
    }
  }
  int possible_peak_num = std::ceil(win_size * 100);

  double n = possible_peak_num;
  double k = env_peak_intensities.size();
  double n1 = inte_peak_num;
  double k1 = match_peak_num;
  double lambda = n1 / n * k;

  boost::math::poisson_distribution<> pd(lambda); 

  double pvalue = 1.0 - boost::math::cdf(pd, k1);
  return pvalue;
}

double compRankSumPValue(double n1, double n2, double r1) {
  double u1 = n1 * n2 + n1 * (n1 + 1) * 0.5 - r1;

  double mean_u = 0.5 * (n1 * n2);
  double log_sig_u = 0.5 * (std::log(n1) + std::log(n2) + std::log(n1 + n2 + 1) - std::log(12));
  double sig_u = std::exp(log_sig_u);
  LOG_DEBUG("log sig u " << log_sig_u << " sig u " << sig_u << " n1 " << n1 << " n2 " << n2 << " r1 " << r1);
  // all peaks are matched
  if (sig_u == 0.0) {
    return 0.0;
  }
  boost::math::normal_distribution<> nd(mean_u, sig_u);
  double p_value = boost::math::cdf(nd, u1);
  p_value = std::min(p_value, 1-p_value);
  return std::abs(p_value);
}

double getRankSumPValue(PeakPtrVec win_peaks, std::vector<double> &env_peak_intensities) {
  int rank_sum = 0;
  int match_peak_num = 0;

  std::sort(win_peaks.begin(), win_peaks.end(), Peak::cmpInteDec);

  for (size_t i = 0; i < env_peak_intensities.size(); i++) {
    double peak_inte = env_peak_intensities[i];
    if (peak_inte > 0.0) {
      int rank = win_peaks.size();
      for (size_t j = 0; j < win_peaks.size(); j++) {
        if (peak_inte >= win_peaks[j]->getIntensity()) {
          rank = j+1;
          break;
        }
      }
      rank_sum += rank;
      match_peak_num++; 
    }
  }
  int peak_num = win_peaks.size();
  double pvalue = compRankSumPValue(peak_num, match_peak_num, rank_sum);
  return pvalue;
}

PeakCluster::PeakCluster(EnvelopePtr theo_env) {
  theo_env_ = theo_env;
  rep_mass_ = theo_env_->getMonoNeutralMass();
  rep_charge_ = theo_env_->getCharge();

  int peak_num = theo_env_->getPeakNum();
  rep_summed_intensities_.resize(peak_num, 0.0);

  clearScores();

  flag_ = 0;
  init_score_ = false;
  smoother_ = std::make_shared<SavitzkyGolay>(9, 2);
}

void PeakCluster::addEnvelopes(FracFeaturePtr feature_ptr, 
                               RealEnvPtrVec envs) {

  int row_num = feature_ptr->getMaxCharge() - feature_ptr->getMinCharge() + 1;
  int col_num = feature_ptr->getMaxMs1Id() - feature_ptr->getMinMs1Id() + 1;
  
  min_charge_ = feature_ptr->getMinCharge();
  max_charge_ = feature_ptr->getMaxCharge();

  min_ms1_id_ = feature_ptr->getMinMs1Id();
  max_ms1_id_ = feature_ptr->getMaxMs1Id();

  scan_begin_ = feature_ptr->getScanBegin();
  scan_end_ = feature_ptr->getScanEnd();
  LOG_DEBUG("add env row " << row_num << " col " << col_num);

  real_envs_.resize(row_num);

  for (int i = 0; i < row_num; i++) {
    real_envs_[i].resize(col_num);
  }

  for (size_t i = 0; i < envs.size(); i++) {
    int row = envs[i]->getCharge() - min_charge_;
    int col = envs[i]->getSpId() - min_ms1_id_;
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
  sum_dist_scores_.resize(2, 1.0);
  sum_corr_scores_.resize(2, 0.0);
  sum_inte_scores_.resize(2, 0.0);
  xic_corr_between_best_charges_.resize(2, 0.0);
}

void PeakCluster::updateScore(PeakPtrVec2D &raw_peaks, bool check_pvalue) {
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

    charge_xic[i].resize(xic_len, 0.0);

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
        summed_intensities[k] += env->getIntensity(k);
      }

      int ms1_id = min_ms1_id_ + j;
      PeakPtrVec all_peaks = raw_peaks[ms1_id];
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
        charge_xic[i][xic_start_idx+j] = env_inte_sum;
      }

      bool level_one_env = true;
      bool level_two_env = true;
      if (check_pvalue) {
        double poisson_pvalue = getPoissonPValue(win_peaks, win_size_, real_intensities);
        double rank_sum_pvalue = getRankSumPValue(win_peaks, real_intensities);
        level_one_env = (rank_sum_pvalue < 0.01 && poisson_pvalue < 0.01);
        //levelTwoEnvelope = (rankSumPValue < 0.05 || poissonPValue < 0.05);
      }
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
      sum_dist_scores_[charge_idx] = std::min(sum_dist_scores_[charge_idx], bc_dist);
      double pc = getPearsonCorr(theo_intensities, summed_intensities);
      sum_corr_scores_[charge_idx] = std::max(sum_corr_scores_[charge_idx], pc);

      if (best_charges_[charge_idx] < 1 || bc_dist < best_charge_dists[charge_idx]) {
        best_charges_[charge_idx] = charge;
        best_charge_dists[charge_idx] = bc_dist;
        if (summed_win_high_inte > 0.0) {
          sum_inte_scores_[charge_idx] = summed_iso_high_inte/summed_win_high_inte;
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
    if (even_best_charge >= 0 && odd_best_charge >= 0) {
      std::vector<double> v1 = smoother_->smooth(charge_xic[even_best_charge]);
      std::vector<double> v2 = smoother_->smooth(charge_xic[odd_best_charge]);
      xic_corr_between_best_charges_[0] = getPearsonCorr(v1, v2); 
      v1 = smoother_->smooth(xic2[even_charge_idx_]);
      v2 = smoother_->smooth(xic2[odd_charge_idx_]);
      xic_corr_between_best_charges_[1] = getPearsonCorr(v1, v2);
    }
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

