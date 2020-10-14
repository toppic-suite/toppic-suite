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

#ifndef TOPPIC_TOPFD_ENV_ENV_PARA_HPP_
#define TOPPIC_TOPFD_ENV_ENV_PARA_HPP_

#include <memory>
#include <vector>
#include <string>
#include <map>

#include "topfd/common/topfd_para.hpp"

namespace toppic {
class EnvPara;
typedef std::shared_ptr<EnvPara> EnvParaPtr;

class EnvPara {
 public:
  EnvPara(){};
  EnvPara(TopfdParaPtr topfd_para_ptr); 
  EnvPara(EnvParaPtr env_para_ptr){
    max_charge_ = env_para_ptr->max_charge_;
    max_mass_ = env_para_ptr->max_mass_;
    window_size_ = env_para_ptr->window_size_;
    estimate_min_inte_ = env_para_ptr->estimate_min_inte_;
    ms_two_sn_ratio_ = env_para_ptr->ms_two_sn_ratio_;
    ms_one_sn_ratio_ = env_para_ptr->ms_one_sn_ratio_;
    min_inte_ = env_para_ptr->min_inte_;
    min_refer_inte_ = env_para_ptr->min_refer_inte_;
    mz_tolerance_ = env_para_ptr->mz_tolerance_;
    min_mass_ = env_para_ptr->min_mass_;
    mass_group_boundary_ = env_para_ptr->mass_group_boundary_;
    percentage_bound_ = env_para_ptr->percentage_bound_;
    max_back_peak_num_ = env_para_ptr->max_back_peak_num_;
    max_forw_peak_num_ = env_para_ptr->max_forw_peak_num_;
    min_match_peak_num_ = env_para_ptr->min_match_peak_num_;
    check_consecutive_peak_num_ = env_para_ptr->check_consecutive_peak_num_;
    relative_consecutive_peak_num_ = env_para_ptr-> relative_consecutive_peak_num_;
    min_consecutive_peak_num_ = env_para_ptr-> min_consecutive_peak_num_;
    do_mz_shift_ = env_para_ptr->do_mz_shift_;
    shift_scale_ = env_para_ptr->shift_scale_;
    do_inte_ratio_ = env_para_ptr->do_inte_ratio_;
    inte_ratio_scale_ = env_para_ptr->inte_ratio_scale_;
    bgn_ratio_ = env_para_ptr->bgn_ratio_;
    end_ratio_ = env_para_ptr->end_ratio_;
    score_error_tolerance_ = env_para_ptr->score_error_tolerance_;
    min_match_env_score_ = env_para_ptr->min_match_env_score_;
    charge_computation_bgn_ = env_para_ptr->charge_computation_bgn_;
    charge_computation_mz_tolerance_ = env_para_ptr->charge_computation_mz_tolerance_;
    rank_peak_distance_ = env_para_ptr->rank_peak_distance_;
    max_similar_mz_env_rank_ = env_para_ptr->max_similar_mz_env_rank_;
    env_num_per_window_ = env_para_ptr->env_num_per_window_;
    do_final_filtering_ = env_para_ptr->do_final_filtering_;
    low_high_dividor_ = env_para_ptr->low_high_dividor_;
    aa_avg_mass_ = env_para_ptr->aa_avg_mass_;
    peak_density_ = env_para_ptr->peak_density_;
    keep_unused_peaks_ = env_para_ptr->keep_unused_peaks_;
    output_multiple_mass_ = env_para_ptr->output_multiple_mass_;
    multiple_min_mass_ = env_para_ptr->multiple_min_mass_;
    multiple_min_charge_ =env_para_ptr->multiple_min_charge_;
    multiple_min_ratio_ = env_para_ptr->multiple_min_ratio_;
    prec_deconv_interval_ = env_para_ptr->prec_deconv_interval_;
    use_env_cnn_ = env_para_ptr->use_env_cnn_;
  }; 

  int getMassGroup(double base_mass);

  void setMinInte(double min_inte, int ms_level);

  int compMinConsPeakNum(int peak_num, int mass_group);

  void setTolerance(double tolerance); 

  double getMzTolerance() {return mz_tolerance_;}

  double getScoreErrorTolerance() {return score_error_tolerance_;}

  double getPercentBound(int mass_group) {return percentage_bound_[mass_group];}

  static int getDefaultMaxCharge() {return 30;}

  static double getDefaultMaxMass() {return 100000;}

  // using input parameters to assign: max_chrg, max_mass 
  int max_charge_ = 30;
  double max_mass_ = 100000;
  // window size 1 m/z
  double window_size_ = 1.0;

  // preprocessing
  // estimate min intensity using thrash method. 
  bool estimate_min_inte_ = true;
  // signal noise ratio 
  double ms_two_sn_ratio_ = 1;
  // ms one signal noise ratio
  double ms_one_sn_ratio_ = 3;
  // minimum peak intensity 
  double min_inte_ = 0;
  // minimum base peak intensity 
  // min_refer_inte_ = min_inte * sn_ratio_ 
  double min_refer_inte_ = 0;


  /*************************************
  // Distribution Envelope factory 
  // initialized before deconvolution 

  // EnvBasePtr is moved to the class EnvBase
  //static EnvBasePtr env_base_ptr_;
  //std::string distr_file_name_ = "/base_data/theo_patt.txt";
  //int distr_entry_num_ = 11000;
  //double distr_mass_interval_ = 10;

  //std::string env_rescore_para_file_name_ = "/base_data/env_rescore_para.txt";

  //std::vector<std::vector<double> > env_rescore_para_;
  **/

  // Envelope detection
  // min_inte and min_ref_inte are used in envelope detection

  // error tolerance for matching real peaks to theoretical peaks 
  double mz_tolerance_ = 0.02;

  // the minimum monoisotopic envelope for an envelope 
  double min_mass_ = 50;

  // Several parameters are related with the mass of envelopes. We classify
  // envelopes into 3 groups based on its base mass. See getMassGroup().
  std::vector<double> mass_group_boundary_ = {min_mass_, min_mass_, 1500, max_mass_};

  // perc bound is used for remove low intensities in theoretical envelopes
  std::vector<double> percentage_bound_ = {0.95, 0.95, 0.85};

  //  maximum number of peaks left and right to the base peak in theoretical envelopes
  int max_back_peak_num_ = 8;
  int max_forw_peak_num_ = 8;

  // Envelope filtering

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
  // if small mz shift is used 
  bool do_mz_shift_ = false;

  // when mz shift is used, the minimum shift is 0.001, shift_fold = 1/0.001
  int shift_scale_ = 1000;
  // if intensity ratio is used 
  bool do_inte_ratio_ = false;

  // when intensity ratio is used, the minimum shift is 0.01, ratio_fold =
  // 1/0.01. We enumerate all possible intensity ratios from bgn_ratio to
  // end_ratio
  int inte_ratio_scale_ = 100;
  double bgn_ratio_ = 0.8;
  double end_ratio_ = 1.2;
  // maximum error in computing m/z accuracy 
  double score_error_tolerance_ = mz_tolerance_; 
  // minimum score for matching envelopes 
  double min_match_env_score_ = 0;

  // 3. fitering using envelopes with charge X, 2X, 3X, ... no parameters here
  // 4. filtering by comparing envelopes with similar charge
  int charge_computation_bgn_ = 15;
  double charge_computation_mz_tolerance_ = 0.002;

  // 5. filtering by comparing envelopes with similar mz and same charge
  // filtering on mz 
  double rank_peak_distance_ = 12;
  int max_similar_mz_env_rank_ = 0;

  // Envelope assigned to 1 m/z intervals
  // number of envelopes per window 
  int env_num_per_window_ = 5;

  // DP algorithm
  //DpParaPtr dp_para_ptr_ = std::make_shared<DpPara>();
  /*  
  // Check double increasing when two envelopes overlap 
  bool check_double_increase_ = true;
  std::vector<std::vector<bool>> coexist_table_;

  // maximum number of envelopes sharing one peak 
  int max_env_num_per_peak_ = 2;
  // used in dpB to specify the number of output envelopes 
  int dp_env_num_ = 300;
  // maximum number of vertices per window 
  int max_env_num_per_vertex_ = 10;
  */

  // envelope final filtering
  // use filtering to  only highest peaks. 
  bool do_final_filtering_ = true;
  double low_high_dividor_ = 1500;
  double aa_avg_mass_ = 120;
  double peak_density_ = 2;

  //  unused peaks
  bool keep_unused_peaks_ = false;

  //  multiple mass 
  bool output_multiple_mass_ = false;
  double multiple_min_mass_ = 5000;
  int multiple_min_charge_ = 20;
  double multiple_min_ratio_ = 0.9;

  // PrecDeconv
  double prec_deconv_interval_ = 3.0;

  // Use EnvCNN
  bool use_env_cnn_ = false;
};

} /* namespace */

#endif 
