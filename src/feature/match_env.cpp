//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include <cstddef>
#include <limits>

#include "feature/match_env.hpp" 

namespace prot {

MatchEnv::MatchEnv(int mass_group, EnvelopePtr theo_env_ptr, 
                   RealEnvPtr real_env_ptr):
    mass_group_(mass_group),
    theo_env_ptr_(theo_env_ptr), 
    real_env_ptr_(real_env_ptr) {
    }

// intensity normalization method 
double calcNormInteScr(double intensity) {
  return (double) std::sqrt(intensity);
}

// Scoring function
// compute matching score_ 
void MatchEnv::compScr(FeatureMngPtr mng_ptr) {
  if (mng_ptr->do_mz_shift_) {
    double best_shift = findBestShift(mng_ptr);
    theo_env_ptr_->changeMz(best_shift);
  }

  if (mng_ptr->do_inte_ratio_) {
    double best_ratio = findBestRatio(mng_ptr);
    theo_env_ptr_->changeIntensity(best_ratio);
  }
  score_ = calcScrWithSftRatio(0, 1, mng_ptr->score_error_tolerance_);
}

// search for best m/z shift 
double MatchEnv::findBestShift(FeatureMngPtr mng_ptr) {
  double best_score = - std::numeric_limits<double>::infinity();
  int best_shift = 0;
  // initialize start shift and end shift based on configuration 
  int bgn_shift = (int)std::round(-mng_ptr->mz_tolerance_ * mng_ptr->shift_scale_);
  int end_shift = (int)std::round(mng_ptr->mz_tolerance_ * mng_ptr->shift_scale_);

  for (int s = bgn_shift; s <= end_shift; s++) {
    double tmp_score = calcScrWithSftRatio((double) s / mng_ptr->shift_scale_,
                                           1.0, mng_ptr->score_error_tolerance_);
    if (tmp_score > best_score) {
      best_score = tmp_score;
      best_shift = s;
    }
  }
  return (double) best_shift / mng_ptr->shift_scale_;
}

// search for best intensity ratio 
double MatchEnv::findBestRatio(FeatureMngPtr mng_ptr) {
  double best_score = - std::numeric_limits<double>::infinity();
  int best_ratio = 0;
  // initialize start ratio and end ratio based on configuration 
  int bgn_ratio = (int)std::round(mng_ptr->bgn_ratio_ * mng_ptr->inte_ratio_scale_);
  int end_ratio = (int)std::round(mng_ptr->end_ratio_ * mng_ptr->inte_ratio_scale_);

  for (int r = bgn_ratio; r <= end_ratio; r++) {
    double tmp_score 
        = calcScrWithSftRatio(0, (double) r/ mng_ptr->inte_ratio_scale_, mng_ptr->score_error_tolerance_);
    if (tmp_score > best_score) {
      best_score = tmp_score;
      best_ratio = r;
    }
  }
  return (double) best_ratio / mng_ptr->inte_ratio_scale_;
}

// Calculating the score_ with shift. 
double MatchEnv::calcScrWithSftRatio(double shift, double ratio, double tolerance) {
  double s = 0;
  for (int i = 0; i < real_env_ptr_->getPeakNum(); i++) {
    // here mz_accu >= 0 and inte_scr >= 0 
    double mz_factor = calcMzFactor(i, shift, tolerance);
    double intensity_factor = calcIntensityFactor(i, ratio);
    double peak_score = mz_factor * intensity_factor
        * calcNormInteScr(theo_env_ptr_->getIntensity(i) * ratio);
    s += peak_score;
  }
  return s;
}


// function of mz accuracy 
double MatchEnv::calcMzFactor(int id_x, double shift, double tolerance) {
  double mz_factor;
    if (real_env_ptr_->isExist(id_x)) {
    double dist = std::abs(theo_env_ptr_->getMz(id_x) + shift - real_env_ptr_->getMz(id_x));
    mz_factor = (tolerance - dist) / tolerance;
    if (mz_factor < 0) {
      mz_factor = 0;
    }
  } else {
    mz_factor = 0;

  }
  return mz_factor;
}

// Calculate the score_ two intensities 
double MatchEnv::calcIntensityFactor(double theo_inte, double real_inte) {
  double ratio = theo_inte / real_inte;
  double intensity_factor;
  // notice the special curve here : 2.0* and sqrt() 
  if (ratio > 1.0) {
    intensity_factor = 1.0 - 2.0 * (ratio - 1.0);
  } else {
    intensity_factor = (double) std::sqrt(ratio);
  }
  if (intensity_factor < 0) {
    intensity_factor = 0;
  }
  return intensity_factor;
}

// function of intensity accuracy 
double MatchEnv::calcIntensityFactor(int id_x, double ratio) {
  double factor;
  if (real_env_ptr_->isExist(id_x)) {
    factor = calcIntensityFactor(theo_env_ptr_->getIntensity(id_x) * ratio,
                                 real_env_ptr_->getIntensity(id_x));
  } else {
    factor = 0;
  }
  return factor;
}


double MatchEnv::calcPeakScr(int id_x, double inte_sum, double tolerance) {
  double mz_factor = calcMzFactor(id_x, 0, tolerance);
  double intensity_factor = calcShareInteAccu(id_x, inte_sum);
  double peak_score = mz_factor * intensity_factor
      * calcNormInteScr(theo_env_ptr_->getIntensity(id_x));
  return peak_score;
}

double MatchEnv::calcShareInteAccu(int id_x, double inte_sum) {
  double intensity_factor;
  double theo_intensity = theo_env_ptr_->getIntensity(id_x);
  if (real_env_ptr_->isExist(id_x)) {
    double real_intensity = real_env_ptr_->getIntensity(id_x);
    double share_ratio = theo_intensity / inte_sum;
    double share_intensity = real_intensity * share_ratio;
    intensity_factor = calcIntensityFactor(theo_intensity, share_intensity);
  } else {
    intensity_factor = 0;
  }
  return intensity_factor;
}

void MsalignWriter::write(std::ofstream &file, MatchEnvPtrVec &envs,
                              MsHeaderPtr header_ptr) {
  file << "BEGIN IONS" << std::endl;
  file << "ID=" << header_ptr->getId() << std::endl;
  file << "SCANS=" << header_ptr->getScansString() << std::endl;
  file << "RETENTION_TIME=" << header_ptr->getRetentionTime() << std::endl;
  if (header_ptr->getActivationPtr() != nullptr) {
    file << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
  }
  if (header_ptr->getMsLevel() > 1) {
    file << "MS_ONE_ID=" << header_ptr->getMsOneId() << std::endl;
    file << "MS_ONE_SCAN=" << header_ptr->getMsOneScan() << std::endl;
    file << "PRECURSOR_MZ=" << header_ptr->getPrecMonoMz() << std::endl;
    file << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;
    file << "PRECURSOR_MASS=" <<  header_ptr->getPrecMonoMass() << std::endl;
    file << "PRECURSOR_INTENSITY=" << header_ptr->getPrecInte() << std::endl;
    if (header_ptr->getFeatureId() >= 0) {
      file << "FEATURE_ID=" << header_ptr->getFeatureId() << std::endl;
      file << "FEATURE_INTENSITY=" << header_ptr->getFeatureInte() << std::endl;
    }
  }

  for (size_t i = 0; i < envs.size(); i++) {
    MatchEnvPtr env = envs[i];
    EnvelopePtr theo_env = env->getTheoEnvPtr();
    RealEnvPtr real_env = env->getRealEnvPtr();
    file << real_env->getMonoMass();
    file << "\t" << theo_env->compIntensitySum();
    file << "\t" << theo_env->getCharge();
    file << std::endl;
  }
  file << "END IONS" << std::endl;
  file << std::endl;
}

}
