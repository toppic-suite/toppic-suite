#ifndef PROT_FEATURE_MNG_HPP_
#define PROT_FEATURE_MNG_HPP_

#include "feature/env_base.hpp"

namespace prot {

class FeatureMng {
 public:
  FeatureMng();

  int getMassGroup(double base_mass);

  void setMinInte(double min_inte);

  int compMinConsPeakNum(int peak_num, int mass_group);

  void setTolerance(double tolerance);

  double getPercentBound(int mass_group) {return percentage_bound_[mass_group];}

  static int getDefaultMaxCharge() {return 30;}

  static double getDefaultMaxMass() {return 100000;}

  // using input parameters to assign: max_chrg, max_mass 
  int max_charge_ = 30;
  double max_mass_ = 100000;
  double window_size_ = 1.0;

  // preprocessing
  // estimate min intensity using thrash method. 
  bool estimate_min_inte_ = true;
  // signal noise ratio 
  double sn_ratio_ = 1;
  // minimum peak intensity 
  double min_inte_ = 0;
  // minimum base peak intensity 
  // min_refer_inte_ = min_inte * sn_ratio_ 
  double min_refer_inte_ = 0;

  // Envelope detection
  // min_inte and min_ref_inte is used in envelope detection

  // Distribution Envelope factory 
  // initialized before deconvolution 
  
   static EnvBasePtr env_base_ptr_;
   std::string distr_file_name_ = "/toppic_resources/base_data/theo_patt.txt";
   int distr_entry_num_ =  11000;
   double distr_mass_interval_ = 10;
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

  // error tolerance for matching real peaks to theoretical peaks 
  double mz_tolerance_ = 0.02;

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

  // Envelope assigned to 1 Da intervals
  // number of envelopes per window 
  int env_num_per_window_ = 5;

  // DP algorithm
  // Check double increasing when two envelopes overlap 
  bool check_double_increase_ = true;
  std::vector<std::vector<bool>> coexist_table_;

  // maximum number of envelopes sharing one peak 
  int max_env_num_per_peak_ = 2;
  // used in dpB to specify the number of output envelopes 
  int dp_env_num_ = 300;
  // maximum number of vertices per window 
  int max_env_num_per_vertex_ = 10;

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
  double prec_deconv_interval_ = 2;
};

typedef std::shared_ptr<FeatureMng> FeatureMngPtr;

} /* namespace */

#endif 
