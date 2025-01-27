//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_TOPFD_ENV_ENV_PARA_HPP_
#define TOPPIC_TOPFD_ENV_ENV_PARA_HPP_

#include <memory>
#include <vector>
#include <limits>

#include "topfd/common/topfd_para.hpp"

namespace toppic {

class EnvPara;
typedef std::shared_ptr<EnvPara> EnvParaPtr;

class EnvPara {
 public:
  //EnvPara() to be removed
  EnvPara(){};
  // set mz tolerance 
  EnvPara(double mz_tolerance);  

  double getMzTolerance(int charge); 

  double getScoreErrorTolerance() {return score_error_tolerance_;}


  //***fixed parameters used in env_detect***
  double getPercentBound(int mass_group) {return percentage_bound_[mass_group];}

  int getMassGroup(double base_mass);

  // the minimum monoisotopic envelope for an envelope 
  double min_mass_ = 50;

  // Several parameters are related with the mass of envelopes. We classify
  // envelopes into 3 groups based on its base mass. See getMassGroup()
  std::vector<double> mass_group_boundary_ = {min_mass_, min_mass_, 1500, std::numeric_limits<double>::max()};

  //  maximum number of peaks left and right to the base peak in theoretical envelopes
  int max_back_peak_num_ = 8;
  int max_forw_peak_num_ = 8;

  // perc bound is used for remove low intensities in theoretical envelopes
  std::vector<double> percentage_bound_ = {0.95, 0.95, 0.85};

  // ***Fixed parameters for envelope filtering ***
  int compMinConsPeakNum(int peak_num, int mass_group);

  // 1. filtering based on real envelop an real envelope is valid if 1) peak
  // number >= 3 2) at most 1 missing peak 3) consecutive peak number >=
  // peak_num - 3, 3
  // minimum peak number in an envelope 
  std::vector<int> min_match_peak_num_ = {1, 2, 3 };
  int max_miss_peak_num_ = 1;
  // check consecutive peaks 
  bool check_consecutive_peak_num_ = true;
  int relative_consecutive_peak_num_ = -3;
  std::vector<int> min_consecutive_peak_num_ = { 1, 2, 3 };

  // 2. filtering using score
  // parameters for computing scores of matching envelopes 
  // Optimize the score using small shifts of theoretical m/z values 
  bool do_mz_shift_ = false;
  // when mz shift is used, the minimum shift is 0.001, shift_fold = 1/0.001
  int shift_scale_ = 1000;

  // Optimize the score using a scale ratio of theoretical peak intensities 
  bool do_inte_ratio_ = false;
  // when intensity ratio is used, the minimum shift is 0.01, ratio_fold =
  // 1/0.01. We enumerate all possible intensity ratios from bgn_ratio to
  // end_ratio
  int inte_ratio_scale_ = 100;
  double bgn_ratio_ = 0.8;
  double end_ratio_ = 1.2;

  // maximum m/z error used in msdeconv score computation of match envelopes.  
  double score_error_tolerance_ = 0.02; 
  // minimum score for matching envelopes 
  double min_match_env_score_ = 0;

  // 3. fitering using envelopes with charge Z, 2 times of Z, 3 times of Z, 
  // no parameters here

  // 4. filtering by comparing envelopes with similar charge.
  // If two envelopes uses the same reference peak and their
  // charge states are Z and Z+1, then one envelope is removed. 
  // This step is applied to only high charge state (Z>= 15) envelopes
  int charge_computation_bgn_ = 15;
  double charge_computation_mz_tolerance_ = 0.002;

  // 5. filtering by comparing envelopes with similar mz and same charge
  // If two envelopes have the same charge and their monoisotopic mass 
  // difference is less than 12 Dalton, then only the top scoring one is keep.
  double rank_peak_distance_ = 12;
  // Only keep the top envelope
  int max_similar_mz_env_rank_ = 0;

  // ***** Fixed parameters for match env filter ***
  // Monoisotopic masses are divided into two groups 
  // < 1500 and > 1500 in filtering
  double low_high_dividor_ = 1500;
  double aa_avg_mass_ = 120;
  double peak_density_ = 2;
  int compLowMassNum();
  int compHighMassNum(double prec_mass);

  // **** Fixed parameters in match_env_util::addMultiMass
  // For one envelope, if we cannot determine its
  // charge state and monoisotopic mass, we will 
  // add several envelopes with two consecutive charges
  // and envelopes with -1 and +1 Dalton shift
  // See match_env_util::addMultipleMass
  double multiple_min_mass_ = 5000;
  int multiple_min_charge_ = 20;
  double multiple_min_ratio_ = 0.9;

 private:
  // *** mz_tolerance and score_error_tolerance,
  // *** are initialized and fixed for envelope detection, filtering, scoring
  // error tolerance for matching real peaks to theoretical peaks 
  double mz_tolerance_ = 0.02;
};

} /* namespace */

#endif 
